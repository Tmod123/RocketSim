import sim2
import matplotlib.pyplot as plt 
import numpy as np
import csv
import random
import json
import os
from multiprocessing import Pool, cpu_count
from matplotlib.patches import Ellipse

N_GRID = 600
 
# (low pct, high pct, label, color, alpha, linewidth) on distance-from-center rank
BUCKETS = [
    (0.00, 0.50, "Inner 50%", "tab:green",  0.15, 0.8),
    (0.50, 0.90, "50-90%",    "tab:orange", 0.25, 0.8),
    (0.90, 1.01, "Outer 10%", "tab:red",    0.60, 1.0),
]
 
 
# Where things sit inside one run's history, and which column is "up".
# Edit these if your history is ordered differently.
T_IDX, POS_IDX, VEL_IDX = 0, 1, 2   # index of time, position, velocity
Z_COL = 2                           # columns are x, y, z
 
 
def _column(a, col):
    """Pull one component from an (N,3) array (or (3,N); or pass through 1-D)."""
    a = np.asarray(a, dtype=float)
    if a.ndim == 0:
        raise ValueError("got a single number where an array was expected")
    if a.ndim == 1:
        return a
    if a.shape[1] != 3 and a.shape[0] == 3:
        a = a.T
    return a[:, col]
 
 
def _describe(run):
    lines = []
    for i, x in enumerate(run):
        try:
            shape = np.asarray(x).shape
        except Exception:
            shape = "?"
        lines.append(f"  [{i}] {type(x).__name__}, shape {shape}")
    return "\n".join(lines)
 
 
def _clean(run):
    """history -> (time, altitude, vertical velocity) as 1-D arrays."""
    try:
        t = np.asarray(run[T_IDX], dtype=float)
        alt = _column(run[POS_IDX], Z_COL)
        vz = _column(run[VEL_IDX], Z_COL)
        if not (len(t) == len(alt) == len(vz)):
            raise ValueError(f"lengths differ: time {len(t)}, position {len(alt)}, velocity {len(vz)}")
    except (IndexError, ValueError, TypeError) as e:
        raise ValueError(
            f"Could not read a history using T_IDX={T_IDX}, POS_IDX={POS_IDX}, VEL_IDX={VEL_IDX} ({e}).\n"
            f"The history has {len(run)} items:\n{_describe(run)}\n"
            "Set T_IDX / POS_IDX / VEL_IDX at the top of this section to match."
        ) from e
    return t, alt, vz
 
 
def _resample(runs, nominal, n_grid=N_GRID):
    """Common time grid; values after a run lands are NaN so its line stops."""
    t_end = max(max(r[0][-1] for r in runs), nominal[0][-1])
    grid = np.linspace(0, t_end, n_grid)
 
    def one(run):
        t, alt, vz = run
        return (np.interp(grid, t, alt, right=np.nan),
                np.interp(grid, t, vz, right=np.nan))
 
    A, V = zip(*[one(r) for r in runs])
    na, nv = one(nominal)
    return grid, np.array(A), np.array(V), na, nv
 
 
def _distance_percentile(runs):
    """Normalized distance of (apogee, max velocity, flight time) from the
    medians, as a percentile rank: 0 = most typical, 1 = most unusual."""
    feats = np.array([[r[1].max(), r[2].max(), r[0][-1]] for r in runs])
    sd = feats.std(axis=0)
    sd[sd == 0] = 1.0
    dist = np.linalg.norm((feats - np.median(feats, axis=0)) / sd, axis=1)
    return dist.argsort().argsort() / max(len(dist) - 1, 1)
 
 
def _panels():
    fig, axes = plt.subplots(1, 2, figsize=(13, 5))
    for ax, ylab in zip(axes, ["Altitude (m)", "Vertical velocity (m/s)"]):
        ax.set_xlabel("Time (s)"); ax.set_ylabel(ylab); ax.grid(alpha=0.3)
    return fig, axes
 
 
def make_ensemble_plots(nominal, runs, outdir=".", show=True):
    """show=True leaves the figures open so your later plt.show() displays them."""
    nominal = _clean(nominal)
    runs = [_clean(r) for r in runs]
    grid, A, V, na, nv = _resample(runs, nominal)
    pct = _distance_percentile(runs)
    data = [(A, na), (V, nv)]
 
    # --- all runs ---
    fig, axes = _panels()
    for ax, (arr, nom) in zip(axes, data):
        for row in arr:
            ax.plot(grid, row, color="tab:blue", alpha=0.08, lw=0.8)
        ax.plot(grid, nom, color="black", lw=2, label="Nominal")
        ax.legend()
    axes[0].set_title(f"Altitude, all runs (n={len(runs)})")
    axes[1].set_title(f"Vertical velocity, all runs (n={len(runs)})")
    fig.tight_layout()
    fig.savefig(os.path.join(outdir, "ensemble_all.png"), dpi=150)
    if not show:
        plt.close(fig)
 
    # --- bucketed ---
    fig, axes = _panels()
    for ax, (arr, nom) in zip(axes, data):
        for lo, hi, label, color, alpha, lw in BUCKETS:   # outliers drawn last
            idx = np.where((pct >= lo) & (pct < hi))[0]
            for k, i in enumerate(idx):
                ax.plot(grid, arr[i], color=color, alpha=alpha, lw=lw,
                        label=f"{label} (n={len(idx)})" if k == 0 else None)
        ax.plot(grid, nom, color="black", lw=2, label="Nominal")
        ax.legend()
    axes[0].set_title("Altitude, colored by distance from center")
    axes[1].set_title("Vertical velocity, colored by distance from center")
    fig.tight_layout()
    fig.savefig(os.path.join(outdir, "ensemble_buckets.png"), dpi=150)
    if not show:
        plt.close(fig)
 
    apo = np.array([r[1].max() for r in runs])
    
 

EAST_COL, NORTH_COL = 0, 1
 
 
def _apogee_and_landing(run):
    pos = np.asarray(run[POS_IDX], dtype=float)
    if pos.shape[1] != 3 and pos.shape[0] == 3:
        pos = pos.T
    k = int(np.argmax(pos[:, Z_COL]))
    ap = (pos[k, EAST_COL], pos[k, NORTH_COL])
    land = (pos[-1, EAST_COL], pos[-1, NORTH_COL])
    return ap, land
 
 
def _draw_ellipses(ax, pts, color, sigmas=(1, 2, 3)):
    """Covariance ellipses around a cloud of (east, north) points."""
    mean = pts.mean(axis=0)
    vals, vecs = np.linalg.eigh(np.cov(pts.T))        # ascending eigenvalues
    angle = np.degrees(np.arctan2(vecs[1, 1], vecs[0, 1]))   # major axis direction
    for n in sigmas:
        w, h = 2 * n * np.sqrt(vals[1]), 2 * n * np.sqrt(vals[0])
        ax.add_patch(Ellipse(mean, w, h, angle=angle, fill=False,
                             edgecolor="black", lw=1.2, zorder=3))
    ax.add_patch(Ellipse(mean, w, h, angle=angle, facecolor=color,
                         alpha=0.25, zorder=2))        # shade the outer (3 sigma)
 
 
def make_dispersion_plot(nominal, runs, outdir=".", show=True,
                         background=None, extent=None, nominal_land=None, nominal_apogee=None):
    """
    background: path to a map/satellite image (optional)
    extent:     (east_min, east_max, north_min, north_max) in meters, relative to
                the launch point, describing what area the image covers.
                Required if background is given.
    nominal_land:   (east, north) of the nominal landing point in meters (optional)
    nominal_apogee:   (east, north) of the nominal apogee point in meters (optional)
    """
    ap = np.array([_apogee_and_landing(r)[0] for r in runs])
    land = np.array([_apogee_and_landing(r)[1] for r in runs])
    pos0 = np.asarray(nominal[POS_IDX], dtype=float)
    launch = (pos0[0, EAST_COL], pos0[0, NORTH_COL]) if pos0.shape[1] == 3 else (0, 0)

    fig, ax = plt.subplots(figsize=(9, 7))
    if background is not None:
        if extent is None:
            raise ValueError("extent=(east_min, east_max, north_min, north_max) is required with background")
        ax.imshow(plt.imread(background), extent=extent, zorder=0, aspect="auto")
 
    ax.scatter(*ap.T, marker="^", s=12, color="green", label="Simulated Apogee", zorder=4)
    ax.scatter(*land.T, marker="v", s=12, color="blue", label="Simulated Landing Point", zorder=4)
    _draw_ellipses(ax, ap, "green")
    _draw_ellipses(ax, land, "blue")
    ax.scatter(*launch, marker="*", s=80, color="black", label="Launch Point", zorder=5)

    ax.scatter(*nominal_apogee, marker="o", s=25, color="orange", label="Nominal Apogee Point", zorder=5)
    ax.scatter(*nominal_land, marker="o", s=25, color="red", label="Nominal Landing Point", zorder=5)

    ax.axhline(launch[1], color="black", lw=0.6, zorder=1)
    ax.axvline(launch[0], color="black", lw=0.6, zorder=1)
    if extent is not None:
        ax.set_xlim(extent[0], extent[1]); ax.set_ylim(extent[2], extent[3])
    else:
        ax.set_aspect("equal", adjustable="datalim"); ax.autoscale_view()
    ax.set_xlabel("East (m)"); ax.set_ylabel("North (m)")
    ax.set_title(r"1$\sigma$, 2$\sigma$ and 3$\sigma$ Dispersion Ellipses: Apogee and Landing Points")
    ax.legend(loc="upper right")
    fig.tight_layout()
    fig.savefig(os.path.join(outdir, "dispersion_ellipses.png"), dpi=150)
    if not show:
        plt.close(fig)
 
    d = np.linalg.norm(land - land.mean(axis=0), axis=1)
    
def import_Motor_Data(filename):
    times = [0.000]
    thrusts = [0.000]
    with open(filename, newline="") as f:
        reader = csv.reader(f)
        for row in reader:
            if not row:
                continue
            try:
                if float(row[0]) != 0.000:
                    times.append(float(row[0]))
                    thrusts.append(float(row[1]))
            except (ValueError, IndexError):
                pass
        return [times, thrusts]


#json importer

def parse_auto_value(v):
    """Handles OpenRocket 'auto X' text values (e.g. bodytube radius)."""
    if isinstance(v, str) and "auto" in v:
        return float(v.split()[-1])
    return float(v)
 
def as_list(x):
    """Normalize a value that may be a single dict or a list of dicts."""
    if x is None:
        return []
    return x if isinstance(x, list) else [x]

def merge_bodytube_subcomponents(bodytubes):
    """A rocket may have one body tube or several (e.g. an upper section
    plus a lower fin/motor section). This merges every body tube's
    subcomponents into one dict of lists, so the rest of the code can
    search across all tubes without caring how many there are."""
    merged = {}
    for bt in bodytubes:
        for key, value in bt.get("subcomponents", {}).items():
            merged.setdefault(key, []).extend(as_list(value))
    return merged

def find_finset(merged_comps):
    """Finds whichever fin-set type is present (trapezoidfinset,
    freeformfinset, ellipticalfinset, ...) across all body tubes.
    Assumes a single fin set overall."""
    for key, values in merged_comps.items():
        if key.endswith("finset"):
            return key, values[0]
    raise KeyError("No fin set found on any body tube.")

def fin_geometry(finset_key, finset):
    """Returns (rootChord, tipChord, sweepLength, finSpan) regardless of
    whether the fin set is trapezoidal (stored directly) or freeform
    (derived from its outline points)."""
    if finset_key == "trapezoidfinset":
        return (
            finset["rootchord"],
            finset["tipchord"],
            finset["sweeplength"],
            finset["height"],
        )
 
    # freeformfinset (and anything else without direct chord/sweep fields):
    # derive an equivalent trapezoid from the fin's outline points.
    points = as_list(finset["finpoints"]["point"])
    xs = [p["@attributes"]["x"] for p in points]
    ys = [p["@attributes"]["y"] for p in points]
 
    span = max(ys)
    root_xs = [x for x, y in zip(xs, ys) if y == min(ys)]
    tip_xs = [x for x, y in zip(xs, ys) if y == span]
 
    root_chord = max(root_xs) - min(root_xs)
    tip_chord = max(tip_xs) - min(tip_xs)  # 0 for a pointed tip
    sweep_length = min(tip_xs) - min(root_xs)
 
    return root_chord, tip_chord, sweep_length, span


def find_launch_guide(merged_comps):
    """Returns (guide_type, length, outer_diameter, inner_diameter).
    Prefers a launch lug (rod-guided); falls back to rail button(s)
    (rail-guided) if no lug is present anywhere. NOTE: a rail button
    doesn't have a true "length" like a lug tube — its height above the
    body is used as the closest analog. These are physically different
    launch systems, so check guideType before assuming rod-launch physics
    apply."""
    if "launchlug" in merged_comps:
        lug = merged_comps["launchlug"][0]
        outer_d = lug["radius"] * 2
        inner_d = (lug["radius"] - lug["thickness"]) * 2
        return "launchlug", lug["length"], outer_d, inner_d
 
    if "railbutton" in merged_comps:
        button = merged_comps["railbutton"][0]  # first button as representative
        return "railbutton", button["height"], button["outerdiameter"], button["innerdiameter"]
 
    raise KeyError("No launchlug or railbutton found on any body tube.")

def find_motor_mount(bodytubes, merged_comps):
    """Motor mount can sit directly on a body tube, or be nested inside an
    inner tube (which itself can live on any body tube). Returns the
    motor mount dict."""
    for bt in bodytubes:
        if "motormount" in bt:
            return bt["motormount"]
 
    for innertube in merged_comps.get("innertube", []):
        if "motormount" in innertube:
            return innertube["motormount"]
 
    raise KeyError("No motor mount found on any body tube or inner tube.")

json_path = "rockets/JSON/rocket.json"

with open(json_path) as f:
    data = json.load(f)
 
rocket = data["rocket"]
stage = rocket["subcomponents"]["stage"]  # single-stage rocket
 
comps = stage["subcomponents"]
nosecone = as_list(comps["nosecone"])[0]  # assumes one nose cone
bodytubes = as_list(comps["bodytube"])    # one tube, or several — handled the same way
 
bt_comps = merge_bodytube_subcomponents(bodytubes)
 
# --- Fins -----------------------------------------------------------
finset_key, finset = find_finset(bt_comps)
rootChord, tipChord, sweepLength, finSpan = fin_geometry(finset_key, finset)
 
numberOfFins = finset["fincount"]
finThickness = finset["thickness"]
jointAngle = finset["cant"]
leadingEdgeAngle = np.arctan(sweepLength / finSpan) if finSpan else 0.0
 
# --- Body -----------------------------------------------------------
noseLength = nosecone["length"]
bodytubeLength = sum(bt["length"] for bt in bodytubes)
# Widest body tube radius — some designs taper across sections
bodyRadius = max(parse_auto_value(bt["radius"]) for bt in bodytubes)
bodyDiameter = bodyRadius * 2
rocketLength = noseLength + bodytubeLength
# Assumes the fin set sits at the very bottom of the last body tube
leadingEdgeToNosecone = rocketLength - rootChord
noseType = nosecone["shape"]
 
# --- Launch guide -----------------------------------------------------
guideType, launchGuideLength, launchGuideOuterDiameter, launchGuideInnerDiameter = (
    find_launch_guide(bt_comps)
)
 
# --- Motor mount -----------------------------------------------------
motormount = find_motor_mount(bodytubes, bt_comps)
motors = as_list(motormount["motor"])
motorDiameter = motors[0]["diameter"]

rodList = []
apogeeList = []
maxVList = []
timeToMaxList = []
totalTimeList = []
groundHitList = []
rangeList = []
machList = []






motor_csv = import_Motor_Data("Motors/Hypertek_L550.csv")
timeCurve = np.array(motor_csv[0])
thrustCurve = np.array(motor_csv[1])











def randomize(value, PercentRange):
    return random.uniform(value*(1-PercentRange), value*(1+PercentRange))

def randomizeList(arr, PercentRange):
    multipliers = np.random.uniform(1 - PercentRange, 1 + PercentRange, size=arr.shape)
    return arr * multipliers

initialTemperature = 20.53583333
TemperaturePercent = .05
initialAltitude = 492
AltitudePercent = .05
initialPressure = 101550
PressurePercent = .02
avgWindSpeed = 5
WindPercent = 0.2
turbulence = 0.15
TurbulencePercent = 0.2
windHeading = np.radians(45)
WindDirectionPercent = 0.15
parachuteDragCoefficient = 0.80
ParaCDPercent = 0.05
parachuteArea = 0.07306
ParaAreaPercent = 0.05
MainDeploymentAltitude = 305
DeploymentPercent = 0.05
MainArea = 1.169
MainAreaPercent = 0.1
MainDragCoefficient = 1.550
MainDragPercent = 0.1
surfaceRoughness = 20e-6
RoughnessPercent = 0.2
launchGuideRoughness = 60e-6
GuideRoughnessPercent = 0.2
dryMass = 4.991
DryPercent = 0.05
fuelMass = 1.552
FuelPercent = 0.05
launchAngle = np.radians(2)
LaunchAnglePercent = 0.01
launchDirection = np.radians(-135)
launchDirectionPercent = 0.10
CG_dry = 1.0463
CGDPercent = 0.02
CG_wet = 1.0774
CGWPercent = 0.02
bodyDiameterPercent = 0.01
lengthPercent = 0.02
FinDimensionsPercent = 0.02
thrustPercent = 0.05
burnPercent = 0.03
sampleCount = 1000


def run_Randoms(simulationNumber, nominal=False):

    # print(
    #     f"Simulation {simulationNumber} | "
    #     f"PID: {os.getpid()} | "
    #     f"Random: {random.random()}"
    # )

    if nominal == True:
        return sim2.run_simulation(initialTemperature, initialAltitude, initialPressure, avgWindSpeed, turbulence, windHeading, parachuteDragCoefficient, parachuteArea, MainDeploymentAltitude, MainArea, MainDragCoefficient, surfaceRoughness, launchGuideRoughness, dryMass, fuelMass, launchAngle, launchDirection, CG_dry, CG_wet, bodyDiameter, thrustCurve, timeCurve, finSpan, rootChord, tipChord, sweepLength, noseLength, bodytubeLength, finThickness, motorDiameter, launchGuideLength, launchGuideOuterDiameter, launchGuideInnerDiameter, numberOfFins)

    initialTemperatureNew = randomize(initialTemperature, TemperaturePercent)
    initialAltitudeNew = randomize(initialAltitude, AltitudePercent)
    initialPressureNew = randomize(initialPressure, PressurePercent)
    avgWindSpeedNew = randomize(avgWindSpeed, WindPercent)
    turbulenceNew = randomize(turbulence, TurbulencePercent)
    windHeadingNew = randomize(windHeading, WindDirectionPercent)
    parachuteDragCoefficientNew = randomize(parachuteDragCoefficient, ParaCDPercent)
    parachuteAreaNew = randomize(parachuteArea, ParaAreaPercent)
    MainDeploymentAltitudeNew = randomize(MainDeploymentAltitude, DeploymentPercent)
    MainAreaNew = randomize(MainArea, MainAreaPercent)
    MainDragCoefficientNew = randomize(MainDragCoefficient, MainDragPercent)
    surfaceRoughnessNew = randomize(surfaceRoughness, RoughnessPercent)
    launchGuideRoughnessNew = randomize(launchGuideRoughness, GuideRoughnessPercent)
    dryMassNew = randomize(dryMass, DryPercent)
    fuelMassNew = randomize(fuelMass, FuelPercent)
    launchAngleNew = randomize(launchAngle, LaunchAnglePercent)
    launchDirectionNew = randomize(launchDirection, launchDirectionPercent)
    CG_dryNew = randomize(CG_dry, CGDPercent)
    CG_wetNew = randomize(CG_wet, CGWPercent)

    bodyDiameterNew = randomize(bodyDiameter, bodyDiameterPercent)
    finSpanNew = randomize(finSpan, FinDimensionsPercent)
    rootChordNew = randomize(rootChord, FinDimensionsPercent)
    tipChordNew = randomize(tipChord, FinDimensionsPercent)
    sweepLengthNew = randomize(sweepLength, FinDimensionsPercent)
    noseLengthNew = randomize(noseLength, lengthPercent)
    bodytubeLengthNew = randomize(bodytubeLength, lengthPercent)
    thrustCurveNew = randomizeList(thrustCurve, thrustPercent)
    timeCurveNew = randomizeList(timeCurve, burnPercent)
    finThicknessNew = randomize(finThickness, FinDimensionsPercent)                                                                                                                                                                                                                                                                                                                                                    


    sort_indices = np.argsort(timeCurveNew)
    timeCurveNew = timeCurveNew[sort_indices]
    thrustCurveNew = thrustCurveNew[sort_indices]
    



    Result =  sim2.run_simulation(initialTemperatureNew, initialAltitudeNew, initialPressureNew, avgWindSpeedNew, turbulenceNew, windHeadingNew, parachuteDragCoefficientNew, parachuteAreaNew, MainDeploymentAltitudeNew, MainAreaNew, MainDragCoefficientNew, surfaceRoughnessNew, launchGuideRoughnessNew, dryMassNew, fuelMassNew, launchAngleNew, launchDirectionNew, CG_dryNew, CG_wetNew, bodyDiameterNew, thrustCurveNew, timeCurveNew, finSpanNew, rootChordNew, tipChordNew, sweepLengthNew, noseLengthNew, bodytubeLengthNew, finThicknessNew, motorDiameter, launchGuideLength, launchGuideOuterDiameter, launchGuideInnerDiameter, numberOfFins)

    return Result

def Histo_results():
    # Graph 1
    plt.figure(figsize=(10, 6))
    plt.hist(rodList, bins='auto', alpha=0.7, label='Rod Velocity')
    plt.xlabel('Rod Velocity (m/s)')
    plt.ylabel('Frequency')
    plt.title('Rod Velocity Distribution')
    plt.grid(True)


    # Graph 2
    plt.figure(figsize=(10, 6))
    plt.hist(apogeeList, bins='auto', alpha=0.7, label='Apogee')
    plt.xlabel('Apogee (m)')
    plt.ylabel('Frequency')
    plt.title('Apogee Distribution')
    plt.grid(True)


    # Graph 3
    plt.figure(figsize=(10, 6))
    plt.hist(maxVList, bins='auto', alpha=0.7, label='Max Velocity')
    plt.xlabel('Max Velocity (m/s)')
    plt.ylabel('Frequency')
    plt.title('Max Velocity Distribution')
    plt.grid(True)


    # Graph 4
    plt.figure(figsize=(10, 6))
    plt.hist(timeToMaxList, bins='auto', alpha=0.7, label='Time to Apogee')
    plt.xlabel('Time to Apogee (s)')
    plt.ylabel('Frequency')
    plt.title('Time to Apogee Distribution')
    plt.grid(True)


    # Graph 5
    plt.figure(figsize=(10, 6))
    plt.hist(totalTimeList, bins='auto', alpha=0.7, label='Total Flight Time')
    plt.xlabel('Total Flight Time (s)')
    plt.ylabel('Frequency')
    plt.title('Total Flight Time Distribution')
    plt.grid(True)


    # Graph 6
    plt.figure(figsize=(10, 6))
    plt.hist(groundHitList, bins='auto', alpha=0.7, label='Ground Hit Velocity')
    plt.xlabel('Ground Hit Velocity (m/s)')
    plt.ylabel('Frequency')
    plt.title('Ground Hit Velocity Distribution')
    plt.grid(True)


    # Graph 7
    plt.figure(figsize=(10, 6))
    plt.hist(rangeList, bins='auto', alpha=0.7, label='Range')
    plt.xlabel('Range (m)')
    plt.ylabel('Frequency')
    plt.title('Range Distribution')
    plt.grid(True)


    # Graph 8
    plt.figure(figsize=(10, 6))
    plt.hist(machList, bins='auto', alpha=0.7, label='Mach')
    plt.xlabel('Mach')
    plt.ylabel('Frequency')
    plt.title('Mach Distribution')
    plt.grid(True)


    plt.show()


if __name__ == "__main__":
    num_cores = max(1, cpu_count() - 1)  # Use all but one core
    print(f"Using {num_cores} CPU cores")
    print(f"Running {sampleCount} simulations...")

    with Pool(processes=num_cores) as pool:
        results = []
        for i, result in enumerate(pool.imap_unordered(run_Randoms, range(sampleCount)), 1):
            results.append(result)
            print(f"Completed {i}/{sampleCount} simulations", end="\r")

    scalars = [r[:8] for r in results]
    histories = [r[8:11] for r in results]
    (rodList, apogeeList, maxVList, timeToMaxList, totalTimeList, groundHitList, rangeList, machList) = map(list, zip(*scalars))
    nominal_history = run_Randoms(None, nominal=True)[8:11]
    nominalTester = run_Randoms(None, nominal=True)

    make_ensemble_plots(nominal_history, histories, outdir=".", show=True)
    make_dispersion_plot(nominal_history, histories, outdir=".", show=True,
                     background="site.png",
                     extent=(-1115, 885, -1015, 985),
                     nominal_land=(nominal_history[1][0][-1], nominal_history[1][1][-1]),
                     nominal_apogee=(nominal_history[1][0][nominal_history[1][2].argmax()], nominal_history[1][1][nominal_history[1][2].argmax()]))  
    Histo_results()
