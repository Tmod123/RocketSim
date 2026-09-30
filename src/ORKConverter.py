#!/usr/bin/env python3
"""
ork_to_json.py — Convert an OpenRocket (.ork) design file to JSON.

An .ork file is XML — either plain XML text, or (more commonly, when saved
by newer OpenRocket versions) a zip archive containing a single XML entry.
This script handles both.

Usage:
    python3 ork_to_json.py input.ork output.json
    python3 ork_to_json.py input.ork output.json --keep-nan   # emit NaN literals
                                                                 instead of null

What you get in the JSON:
- The full rocket design tree (stages, nose cones, body tubes, fins, motor
  mounts, parachutes, materials, appearance, etc.) exactly as nested in the
  file.
- Every <simulation>, with its launch conditions and <flightdata>.
- Each <databranch>'s flight-data samples, EXPANDED: the raw
  "<datapoint>1.0,2.3,...</datapoint>" comma rows are turned into a real
  JSON object per sample, e.g. {"Time": 1.0, "Altitude": 2.3, ...}, using
  the column names in that databranch's `types="..."` attribute. The
  original column order is also kept as "_columns" on the databranch so
  you don't lose that information.
- All values are type-converted: numbers become JSON numbers, true/false
  become JSON booleans, and NaN becomes JSON null by default (pass
  --keep-nan to instead emit the non-standard NaN literal, which Python's
  json module can read back but which is NOT valid per the JSON spec).

Design notes on the generic XML->JSON mapping:
- Element attributes go under a "@attributes" key.
- Element text (for leaf elements with no children) becomes the value
  directly, type-converted.
- If an element has attributes AND text, the text goes under "#text".
- Repeated child tags (e.g. multiple <stage> under <subcomponents>) become
  a JSON list; a single occurrence stays a plain object/value under its
  tag name. This mirrors the way OpenRocket itself treats these elements.
"""

import sys
import json
import zipfile
import argparse
import xml.etree.ElementTree as ET


# ---------------------------------------------------------------------------
# Value conversion helpers
# ---------------------------------------------------------------------------

def convert_scalar(text, keep_nan=False):
    """Convert a raw XML text value to a Python bool/int/float/str/None."""
    if text is None:
        return None
    t = text.strip()
    if t == "":
        return None
    low = t.lower()
    if low == "true":
        return True
    if low == "false":
        return False
    if low == "nan":
        return float("nan") if keep_nan else None
    if low in ("inf", "+inf", "infinity"):
        return float("inf") if keep_nan else None
    if low in ("-inf", "-infinity"):
        return float("-inf") if keep_nan else None
    # Try int first (avoids turning "5" into 5.0)
    try:
        return int(t)
    except ValueError:
        pass
    try:
        return float(t)
    except ValueError:
        pass
    return t


# ---------------------------------------------------------------------------
# Special handling for <databranch>/<datapoint> flight-data rows
# ---------------------------------------------------------------------------

def parse_databranch(elem, keep_nan=False):
    """
    Convert a <databranch types="Time,Altitude,..."> element into:
        {
          "@attributes": {...},   # name, optimumAltitude, etc. (minus 'types')
          "_columns": ["Time", "Altitude", ...],
          "events": [ {...}, ... ],
          "datapoints": [ {"Time": 0.0, "Altitude": 0.0, ...}, ... ]
        }
    """
    attrs = dict(elem.attrib)
    types_str = attrs.pop("types", "")
    columns = [c.strip() for c in types_str.split(",")] if types_str else []

    result = {}
    if attrs:
        result["@attributes"] = {k: convert_scalar(v, keep_nan) for k, v in attrs.items()}
    result["_columns"] = columns

    events = []
    datapoints = []
    for child in elem:
        if child.tag == "event":
            events.append(elem_to_obj(child, keep_nan))
        elif child.tag == "datapoint":
            raw = (child.text or "").split(",")
            if columns and len(raw) == len(columns):
                row = {
                    columns[i]: convert_scalar(raw[i], keep_nan)
                    for i in range(len(columns))
                }
            else:
                # Column count mismatch (shouldn't normally happen) — fall
                # back to a plain list of values so no data is lost.
                row = [convert_scalar(v, keep_nan) for v in raw]
            datapoints.append(row)
        else:
            # Unexpected child type under databranch — keep it generically.
            result.setdefault("_other", []).append(elem_to_obj(child, keep_nan))

    if events:
        result["events"] = events
    result["datapoints"] = datapoints
    return result


# ---------------------------------------------------------------------------
# Generic recursive XML -> JSON-able object conversion
# ---------------------------------------------------------------------------

def elem_to_obj(elem, keep_nan=False):
    if elem.tag == "databranch":
        return parse_databranch(elem, keep_nan)

    children = list(elem)
    obj = {}

    if elem.attrib:
        obj["@attributes"] = {k: convert_scalar(v, keep_nan) for k, v in elem.attrib.items()}

    if children:
        child_map = {}
        for child in children:
            value = elem_to_obj(child, keep_nan)
            tag = child.tag
            if tag in child_map:
                if not isinstance(child_map[tag], list):
                    child_map[tag] = [child_map[tag]]
                child_map[tag].append(value)
            else:
                child_map[tag] = value
        obj.update(child_map)
        return obj

    # Leaf element (no children)
    text_val = convert_scalar(elem.text, keep_nan)
    if obj:  # had attributes -> keep text separately
        if text_val is not None:
            obj["#text"] = text_val
        return obj
    else:
        return text_val


# ---------------------------------------------------------------------------
# .ork file loading (plain XML or zip-wrapped XML)
# ---------------------------------------------------------------------------

def load_ork_root(path):
    """Return the ElementTree root Element of an .ork file, whether it's
    stored as raw XML text or zipped."""
    if zipfile.is_zipfile(path):
        with zipfile.ZipFile(path) as zf:
            xml_names = [n for n in zf.namelist() if n.lower().endswith((".ork", ".xml", ".rkt"))]
            if not xml_names:
                # Fall back to the first file in the archive
                xml_names = zf.namelist()
            with zf.open(xml_names[0]) as f:
                data = f.read()
        return ET.fromstring(data)
    else:
        tree = ET.parse(path)
        return tree.getroot()


def convert_ork_to_json(input_path, output_path, keep_nan=False, indent=2, rocket_only=True):
    keep_nan=False
    indent=2
    rocket_only=True
    root = load_ork_root(input_path)

    if rocket_only:
        # Only convert the <rocket> element (design/geometry/materials).
        # This skips <simulations>, <photostudio>, and <docprefs> entirely.
        rocket_elem = root.find("rocket")
        if rocket_elem is None:
            raise ValueError("No <rocket> element found in the .ork file.")
        data = {"rocket": elem_to_obj(rocket_elem, keep_nan)}
    else:
        data = {root.tag: elem_to_obj(root, keep_nan)}

    with open(output_path, "w") as f:
        json.dump(data, f, indent=indent, allow_nan=keep_nan)
    return data


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main():  
    convert_ork_to_json(
        "Rockets/ORK/rocket.ork",
        "Rockets/JSON/rocket.json"
    )

if __name__ == "__main__":
    main()
