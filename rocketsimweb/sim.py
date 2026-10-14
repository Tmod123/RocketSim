#sim.py - V3.5 (web)
#
# Class-based version of src/sim.py V3.5 for the Streamlit frontend.
# Everything that was a module-level variable in src/sim.py is now an
# attribute of a RocketSim instance (same names, same default values), and
# every function that read those variables is now a method (same names).
# Each browser session gets its own RocketSim, so users never share state.
#
# Pure helper functions that use no simulation state (JSON parsing,
# quaternion math) stay at module level exactly as in src/sim.py.
#
# The rocket design is read from an OpenRocket .ork file through
# ORKConverter (in memory, no JSON file written). A JSON export can still be
# used instead by passing json_path. Relative paths are looked up in the
# folder you run from first, then in this file's folder, so the default
# files are found no matter where Streamlit is started.
#
# Additions for the web app only (none change the physics or the defaults):
#   radius / windSpeed : aliases for bodyRadius / avgWindSpeed, the names
#                        main.py's sidebar sets
#   set_motor_data     : load a thrust curve from an uploaded file
#   windSeed           : optional, makes the wind gusts (and CP jitter) repeatable
#   get_results, results_csv_bytes, get_summary : read the stored histories

import numpy as np
import csv
import io
import os
import random
import json

from ORKConverter import ORKConverter

_HERE = os.path.dirname(os.path.abspath(__file__))


def find_file(path):
    """Use `path` as given if it exists; otherwise look for it relative to
    this file's folder (where the default rocket and motor files live)."""
    if not isinstance(path, str) or os.path.isabs(path) or os.path.exists(path):
        return path
    candidate = os.path.join(_HERE, path)
    return candidate if os.path.exists(candidate) else path


def import_Motor_Data(filename):
    """Same as src/sim.py. `filename` may also be an open text file
    (e.g. io.StringIO of an uploaded file)."""
    times = [0.000]
    thrusts = [0.000]
    if isinstance(filename, str):
        f = open(find_file(filename), newline="")
    else:
        f = filename
        f.seek(0)
    with f:
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


def quaternion_multiply(q1, q2):

    w1, x1, y1, z1 = q1
    w2, x2, y2, z2 = q2

    return np.array([
        w1*w2 - x1*x2 - y1*y2 - z1*z2,
        w1*x2 + x1*w2 + y1*z2 - z1*y2,
        w1*y2 - x1*z2 + y1*w2 + z1*x2,
        w1*z2 + x1*y2 - y1*x2 + z1*w2
    ])


def quaternion_conjugate(q):

    return np.array([
        q[0],
        -q[1],
        -q[2],
        -q[3]
    ])


def quaternion_normalize(q):

    return q / np.linalg.norm(q)


def rotate_vector(q, vector):

    vector_quaternion = np.array([
        0.0,
        vector[0],
        vector[1],
        vector[2]
    ])

    q_conjugate = quaternion_conjugate(q)

    rotated = quaternion_multiply(quaternion_multiply(q, vector_quaternion),q_conjugate)

    return rotated[1:]


def inverse_rotate_vector(q, vector):
    return rotate_vector(quaternion_conjugate(q),vector)


def quaternion_to_euler(q):
    w,x,y,z = q

    pitch = np.arctan2(
        2 * (w*x + y*z),
        1 - 2 * (x*x + y*y)
    )

    yaw = np.arctan2(
        2 * (w*y - z*x),
        1 - 2 * (y*y + x*x)
    )

    roll = np.arctan2(
        2 * (w*z + x*y),
        1 - 2 * (z*z + y*y)
    )

    return np.degrees([pitch,yaw,roll])


def quaternion_from_two_vectors(v1, v2):
    v1 = np.array(v1, dtype=float)
    v2 = np.array(v2, dtype=float)

    v1 = v1 / np.linalg.norm(v1)
    v2 = v2 / np.linalg.norm(v2)

    dot = np.dot(v1, v2)
    cross = np.cross(v1, v2)

    # Vectors point in the same direction
    if dot > 0.999999:
        return np.array([1.0, 0.0, 0.0, 0.0])

    # Vectors point in opposite directions
    if dot < -0.999999:
        axis = np.cross(v1, np.array([1.0, 0.0, 0.0]))

        if np.linalg.norm(axis) < 0.000001:
            axis = np.cross(v1, np.array([0.0, 1.0, 0.0]))

        axis = axis / np.linalg.norm(axis)

        return np.array([
            0.0,
            axis[0],
            axis[1],
            axis[2]
        ])

    quaternion = np.array([
        1.0 + dot,
        cross[0],
        cross[1],
        cross[2]
    ])

    return quaternion_normalize(quaternion)


class RocketSim:

    def __init__(self, ork_path="Rockets/Ork/rocket.ork",
                 motor_path="Motors/Hypertek_L550.csv", json_path=None):
        """ork_path: OpenRocket design (path, bytes, or open binary file).
        json_path: optional - use an existing JSON export instead of ork_path."""
        self.ork_path = ork_path
        self.json_path = json_path

        if json_path is not None:
            with open(find_file(json_path)) as f:
                data = json.load(f)
        else:
            source = find_file(ork_path) if isinstance(ork_path, str) else ork_path
            data = ORKConverter().convert(source)

        rocket = data["rocket"]
        stage = rocket["subcomponents"]["stage"]  # single-stage rocket

        comps = stage["subcomponents"]
        nosecone = as_list(comps["nosecone"])[0]  # assumes one nose cone
        bodytubes = as_list(comps["bodytube"])    # one tube, or several — handled the same way

        bt_comps = merge_bodytube_subcomponents(bodytubes)

        # --- Fins -----------------------------------------------------------
        finset_key, finset = find_finset(bt_comps)
        self.rootChord, self.tipChord, self.sweepLength, self.finSpan = fin_geometry(finset_key, finset)

        self.numberOfFins = finset["fincount"]
        self.finThickness = finset["thickness"]
        self.jointAngle = finset["cant"]
        self.leadingEdgeAngle = np.arctan(self.sweepLength / self.finSpan) if self.finSpan else 0.0

        # --- Body -----------------------------------------------------------
        self.noseLength = nosecone["length"]
        self.bodytubeLength = sum(bt["length"] for bt in bodytubes)
        # Widest body tube radius — some designs taper across sections
        self.bodyRadius = max(parse_auto_value(bt["radius"]) for bt in bodytubes)
        self.bodyDiameter = self.bodyRadius * 2
        self.rocketLength = self.noseLength + self.bodytubeLength
        # Assumes the fin set sits at the very bottom of the last body tube
        self.leadingEdgeToNosecone = self.rocketLength - self.rootChord
        self.noseType = nosecone["shape"]

        # --- Launch guide -----------------------------------------------------
        self.guideType, self.launchGuideLength, self.launchGuideOuterDiameter, self.launchGuideInnerDiameter = (
            find_launch_guide(bt_comps)
        )

        # --- Motor mount -----------------------------------------------------
        motormount = find_motor_mount(bodytubes, bt_comps)
        motors = as_list(motormount["motor"])
        self.motorDiameter = motors[0]["diameter"]

        # --- Parachute(s) -----------------------------------------------------
        # A rocket may have more than one (e.g. drogue + main), possibly on
        # different body tubes. Compute the area for each; the largest-diameter
        # one is treated as "the" main parachute for a single parachuteArea value.
        parachutes = bt_comps.get("parachute", [])
        self.parachuteAreas_all = [np.pi * (p["diameter"] / 2) ** 2 for p in parachutes]
        if self.parachuteAreas_all:
            main_index = max(range(len(parachutes)), key=lambda i: parachutes[i]["diameter"])
            self.parachuteArea = self.parachuteAreas_all[main_index]
        else:
            self.parachuteArea = None
        self.parachuteDragCoefficient = 0.80  # not stored numerically ("auto" in file)
        self.parachuteArea = 0.07306
        self.finSets = 1
        self.transitions = 1 if "transition" in comps else 0
        self.MainDeploymentAltitude = 305
        self.MainArea = 1.169
        self.MainDragCoefficient = 1.550
        self.epsilonAngle = 0 #basically transition angle for nosecone, 0 since smooth, but a cone could be like 15
        self.surfaceRoughness = 20e-6    # meters — needs a materials-roughness table/guess
        self.launchGuideRoughness = 60e-6    # meters — same as above
        self.dryMass = 4.991  # kg — mass properties computed during simulation
        self.fuelMass = 1.552   # kg — motor propellant mass
        self.rodLength   = 1.0                  # m — launch condition, not part of <rocket>
        self.launchAngle = np.radians(2)
        self.launchDirection = np.radians(-135) #heading
        self.launch_vector = np.array([np.sin(self.launchAngle)*np.cos(self.launchDirection), np.sin(self.launchAngle) * np.sin(self.launchDirection), np.cos(self.launchAngle)])
        self.CG_dry = 1.0463  # m from nose tip — computed mass property
        self.CG_wet = 1.0774  # m from nose tip — same

        #Constants
        self.gravity = 9.806 # m/s^2
        self.timeStep = 0.05 # seconds
        self.initialTemperature = 20.53583333 #degrees Celsius
        self.initialAltitude = 492 #meters above sea level
        self.initialPressure = 101550 #Pa
        self.airMolarMass = 0.0289644 #kg/mol
        self.universeGasConstant = 8.31447 #J/(mol*K)
        self.SpecificGasConstant = 287.05 #J/(kg*K)
        self.temperatureLapseRate = 0.0065 #K/m
        self.AdiabaticIndex = 1.4 #unitless & constant for our purpose, value for air
        self.x_n1 = 0
        self.x_n2 = 0
        self.avgWindSpeed = 5 #m/s
        self.turbulence = 0.15 #Turbulence Intensity.
        self.windHeading = np.radians(45) #radians, 0 is wind to the north, pi/2 is wind to the east, pi is wind to the south, 3pi/2 is wind to the west
        self.windSeed = None  # web only: int for repeatable gusts, None = new gusts every run

        # --- CP uncertainty (calibers) ---
        self.cpBiasZ = random.gauss(0,1)   # per-instance draw ~N(0,1); 0 = nominal
        self.cpSigmaBase = 0.10     # cal, 1-sigma away from Mach 1
        self.cpSigmaPeak = 0.30     # cal, extra 1-sigma at Mach 1
        self.cpWidth = 0.22         # Mach half-width of the bump
        self.cpJitter = 0.1         # cal, time-varying wander (0 = off)
        self.cpNoiseTimes, self.cpNoise = [], []

        #motor_csv = import_Motor_Data(fileName)
        self.set_motor_data(motor_path)

        self.calculate_values()
        self.reset_state()

    # ------------------------------------------------------------------
    # Names main.py's sidebar uses, mapped onto the src/sim.py variables
    # ------------------------------------------------------------------
    @property
    def radius(self):
        return self.bodyRadius

    @radius.setter
    def radius(self, value):
        self.bodyRadius = value

    @property
    def windSpeed(self):
        return self.avgWindSpeed

    @windSpeed.setter
    def windSpeed(self, value):
        self.avgWindSpeed = value

    # ------------------------------------------------------------------
    # Motor data
    # ------------------------------------------------------------------
    def set_motor_data(self, source):
        """Load a thrust curve from a path or an open text file (e.g.
        io.StringIO of a Streamlit upload), as src/sim.py does at startup."""
        self.motor_csv = import_Motor_Data(source)
        self.timeCurve = np.array(self.motor_csv[0])
        self.thrustCurve = np.array(self.motor_csv[1])

        self._segment_impulse = 0.5 * (self.thrustCurve[1:] + self.thrustCurve[:-1]) * np.diff(self.timeCurve)
        self.cumulativeImpulse = np.concatenate(([0.0], np.cumsum(self._segment_impulse)))
        self.netImpulse = self.cumulativeImpulse[-1]

    # ------------------------------------------------------------------
    # Values derived from the parameters. src/sim.py computes these once at
    # import; here they are recomputed at the start of every run so edits
    # from the sidebar (radius, rocket length, launch angle, CG, ...) apply.
    # ------------------------------------------------------------------
    def calculate_values(self):
        self.bodyDiameter = self.bodyRadius * 2
        # Assumes the fin set sits at the very bottom of the last body tube
        self.leadingEdgeToNosecone = self.rocketLength - self.rootChord
        self.leadingEdgeAngle = np.arctan(self.sweepLength / self.finSpan) if self.finSpan else 0.0
        self.launch_vector = np.array([np.sin(self.launchAngle)*np.cos(self.launchDirection), np.sin(self.launchAngle) * np.sin(self.launchDirection), np.cos(self.launchAngle)])

        finSpan, sweepLength, tipChord, rootChord = self.finSpan, self.sweepLength, self.tipChord, self.rootChord
        bodyRadius, bodyDiameter, noseLength = self.bodyRadius, self.bodyDiameter, self.noseLength
        rocketLength, numberOfFins, finThickness = self.rocketLength, self.numberOfFins, self.finThickness

        #calculated values
        self.Y1 = np.sqrt(finSpan**2+(sweepLength+tipChord/2-rootChord/2)**2) #supposed to be lambda or the distance on Fin at the mid-chord lines for set 1 (m)
        self.A_fin = (finSpan/2)*(tipChord+rootChord)
        self.finenessRatio = rocketLength/bodyDiameter #unitless, ratio of the length of the rocket to the diameter of the rocket
        self.meanAerodynamicChordLengthOfFin = (2/3) * (rootChord + tipChord - (rootChord*tipChord)/(rootChord + tipChord)) #meters, this is a typical value for a fin
        self.A_wet_nose = (np.pi/(2*bodyRadius**2))*(bodyRadius*noseLength*(bodyRadius**2-noseLength**2)+((bodyRadius**2+noseLength**2)**2)*np.arctan(bodyRadius/noseLength)) #meters^2, this is the wetted area of an ogive nosecone of the rocket
        self.noseVolume = (np.pi/(24*bodyRadius**3))*(6*noseLength*bodyRadius**5+6*bodyRadius*noseLength**5+4*(bodyRadius**3)*(noseLength**3)+6*(bodyRadius**2-noseLength**2)*(bodyRadius**2+noseLength**2)**2*np.arctan(bodyRadius/noseLength))
        self.A_wet_body = np.pi * bodyDiameter * (rocketLength-noseLength) + self.A_wet_nose #meters^2, this is the wetted area of the body of the rocket
        self.A_wet_fins = numberOfFins*(finSpan*(rootChord+tipChord)+finThickness*(rootChord+tipChord+np.sqrt(sweepLength**2 + finSpan**2)+np.sqrt((sweepLength+tipChord-rootChord)**2+finSpan**2))) #meters^2, this is the wetted area of the fins of the rocket
        self.A_ref = np.pi * (bodyRadius**2) #meters^2, this is the reference area of the rocket, which is the cross-sectional area of the rocket body

        self.CP = self.initcenterOfPressure() # m from the nose tip
        self.CT = (self.CG_wet*(self.dryMass+self.fuelMass)-self.CG_dry*self.dryMass)/(self.fuelMass)

    def generate_wind(self):
        rng = random.Random(self.windSeed)  # per-instance, so sessions never share gusts
        self.x_n1 = 0
        self.x_n2 = 0
        self.timeWindSpeeds = []
        self.windSpeeds = []

        for i in range(0,10000):
            w_n = rng.gauss(0,1)
            x_n = w_n + (5/6)*(self.x_n1)-(5/24)*(self.x_n2)
            self.x_n2 = self.x_n1
            self.x_n1 = x_n
            stdev = self.avgWindSpeed*self.turbulence
            self.windSpeeds.append(self.avgWindSpeed + stdev*x_n)
            self.timeWindSpeeds.append(i*0.1)

    def generate_cp_noise(self):
        # Same AR(2) filter as the wind. Seeded off windSeed (+1 so it isn't
        # the identical sequence as the gusts) when a seed is set.
        rng = random.Random(None if self.windSeed is None else self.windSeed + 1)
        self.cpNoiseTimes = [i*0.1 for i in range(10000)]
        self.cpNoise = []
        c1, c2 = 0.0, 0.0
        for i in range(10000):
            c = rng.gauss(0,1) + (5/6)*c1 - (5/24)*c2
            c2, c1 = c1, c
            self.cpNoise.append(c)

    def reset_state(self):
        self.position = np.array([0.0, 0.0, 0.0]) #North, East, Up, is world orientated
        self.velocity = np.array([0.0, 0.0, 0.0]) #North, East, Up, is world orientated
        self.orientation = quaternion_from_two_vectors(np.array([0.0,0.0,1.0]),self.launch_vector) #w, x, y, z. where x is pitch, y is yaw, z is roll, is body orientated
        self.angular_velocity = np.array([0.0, 0.0, 0.0]) #x, y, z, is body orientated
        self.planarHistory = []
        self.positionHistory = np.empty((3,0))
        self.velocityHistory = np.empty((3,0))
        self.forceHistory = np.empty((3,0))
        self.orientationHistory = np.empty((3,0))
        self.angularVelocityHistory = np.empty((3,0))
        self.momentHistory = np.empty((3,0))

        #Init Variables
        self.time = 0
        self.leftRod = False

        #Init Lists
        self.timeList = []
        self.massList = []
        self.machList = []
        self.AOAList = []
        self.CDList = []

    def initcenterOfPressure(self):
        CN_n = 2
        CN_f1 = (1+(self.bodyRadius)/(self.finSpan+self.bodyRadius)) * ((4*self.numberOfFins*(self.finSpan/self.bodyDiameter)**2)/(1+np.sqrt(1+((2*self.Y1)/(self.rootChord+self.tipChord))**2)))
        Pf1 = CN_f1 * (self.leadingEdgeToNosecone+(self.sweepLength*(self.rootChord+2*self.tipChord)/(3*(self.rootChord+self.tipChord)))+(1/6)*(self.rootChord+self.tipChord-(self.rootChord*self.tipChord)/(self.rootChord+self.tipChord)))

        if self.noseType == "cone":
            Pn = CN_n * (0.6667*self.noseLength)
        elif self.noseType == "ogive":
            Pn = CN_n * (0.466*self.noseLength)
        elif self.noseType == "paraboloid":
            Pn = CN_n * (0.5*self.noseLength)
        elif self.noseType == "ellipsoid":
            Pn = CN_n * (0.3333*self.noseLength)
        else:
            Pn = 0

        netNormalForce = CN_n + CN_f1
        netMoment = Pn + Pf1
        return netMoment/netNormalForce

    # ------------------------------------------------------------------
    # Atmosphere
    # ------------------------------------------------------------------
    def temperature(self, position):
        altitude = position[2]
        if altitude+self.initialAltitude < 11000:
            return self.initialTemperature - self.temperatureLapseRate * (altitude)
        else:
            return self.initialTemperature - self.temperatureLapseRate * (11000-self.initialAltitude)

    def pressure(self, position):
        return self.initialPressure * ((self.temperature(position)+273.15)/(self.initialTemperature+273.15))**(self.gravity*self.airMolarMass/(self.universeGasConstant*self.temperatureLapseRate))

    def airDensity(self, position):
        return self.pressure(position)/(self.SpecificGasConstant*(self.temperature(position)+273.15))

    def speedOfSound(self, position):
        return np.sqrt(self.AdiabaticIndex*self.SpecificGasConstant*(self.temperature(position)+273.15))

    def machNumber(self, airVelocity, position):
        return self.speed(airVelocity)/self.speedOfSound(position)

    @staticmethod
    def speed(velocity):
        return np.linalg.norm(velocity)

    @staticmethod
    def angle_of_attack(velocity):
        if velocity[0] == 0 and velocity[1] == 0:
            return 0.0
        else:
            return np.arctan2(np.linalg.norm(velocity[:2]), velocity[2])

    def dynamicViscosity(self, position):
        return 1.458*10**(-6)*((self.temperature(position)+273.15)**(3/2))/(self.temperature(position)+273.15+110.4)

    def kinematicViscosity(self, position):
        return self.dynamicViscosity(position)/self.airDensity(position)

    def reynoldsNumber(self, airVelocity, position, characteristicLength):
        return self.speed(airVelocity) * characteristicLength / self.kinematicViscosity(position)

    # ------------------------------------------------------------------
    # Drag
    # ------------------------------------------------------------------
    def zeroAngleDragCoefficient(self, time, airVelocity, position):
        C_d = self.C_D_friction(airVelocity, position) + self.nosePressureDrag(airVelocity, position) + self.finPressureDrag(time, airVelocity, position) + self.baseDragCoefficient(time, airVelocity, position) + self.parasiticDrag(airVelocity, position)
        return C_d

    def C_D_friction(self, airVelocity, position):
        R = self.reynoldsNumber(airVelocity, position, self.rocketLength)
        M = self.machNumber(airVelocity, position)
        R_crit = 51*(self.surfaceRoughness/self.rocketLength)**-1.039
        if R < 10**4:
            C_f = 1.48*10**-2
        elif R < R_crit:
            C_f = 1/(1.5*np.log(R)-5.6)**2
        else:
            C_f = 0.032*(self.surfaceRoughness/self.rocketLength)**0.2

        C_Mach = (1-0.1*M**2)
        C_f_rough = 0.032*(self.surfaceRoughness/self.rocketLength)**0.2 * C_Mach
        C_f_component = max(C_f, C_f_rough)
        K_body = 1 + 1/(2*self.finenessRatio)
        K_fin = 1 + 2*self.finThickness/self.meanAerodynamicChordLengthOfFin
        C_D_body = C_f_component * K_body * (self.A_wet_body/self.A_ref)
        C_D_fins = C_f_component * K_fin * (self.A_wet_fins/self.A_ref)
        return C_D_body + C_D_fins

    def nosePressureDrag(self, airVelocity, position):
        jointAngle_phi=0.0
        kappa=1.0
        M = self.machNumber(airVelocity, position)
        gamma = 1.4
        eps = np.arctan(self.bodyDiameter/(2*self.noseLength))
        s = np.sin(eps)
        K = 0.72*(kappa-0.5)**2 + 0.82
        C0 = 0.8*np.sin(jointAngle_phi)**2
        C1 = K*s
        dC1 = K*4/(gamma+1)*(1-0.5*s)
        def sup(m):  return K*(2.1*s**2 + 0.5*s/np.sqrt(m**2-1))
        def dsup(m): return -K*0.5*s*m/(m**2-1)**1.5
        if M <= 1:
            a = C1 - C0
            return a*M**(dC1/a) + C0
        if M < 1.3:
            h = 0.3; t = (M-1)/h
            return ((2*t**3-3*t**2+1)*C1 + (t**3-2*t**2+t)*dC1*h
                    + (-2*t**3+3*t**2)*sup(1.3) + (t**3-t**2)*dsup(1.3)*h)
        return sup(M)

    def finPressureDrag(self, time, airVelocity, position):
        M = self.machNumber(airVelocity, position)
        if M < 0.9: #rounded edge
            perpindicularLeadingEdgeDrag = (1-M**2)**(-0.417) -1
        elif M < 1:
            perpindicularLeadingEdgeDrag = 1 - 1.785*(M-0.9)
        else:
            perpindicularLeadingEdgeDrag = 1.214 - (0.502/M**2) + (0.1095/M**4)
        leadingEdgeDrag = perpindicularLeadingEdgeDrag * np.cos(self.leadingEdgeAngle)**2 #angled fin, sweep angle
        trailingEdgeDrag = (1/2) * self.baseDragCoefficient(time, airVelocity, position) #rounded edge
        C_D_fin = leadingEdgeDrag + trailingEdgeDrag
        A_FIN = self.numberOfFins * self.finThickness * self.finSpan
        return (A_FIN/self.A_ref) * C_D_fin

    def baseDragCoefficient(self, time, airVelocity, position):
        if time < self.timeCurve[-1]:
            A_motor = (np.pi/4) * self.motorDiameter ** 2
        else:
            A_motor = 0
        M = self.machNumber(airVelocity, position)
        if M < 1:
            CD_base = 0.12 + 0.13*M**2
        else:
            CD_base = 0.25/M
        return CD_base * (self.A_ref-A_motor)/self.A_ref

    def parasiticDrag(self, airVelocity, position):
        M = self.machNumber(airVelocity, position)
        if M < 1:
            qRatio = 1 + (M**2/4) + (M**4/40)
            C_D_base = 0.12 + 0.13*M**2
        else:
            qRatio = 1.84 - (0.76/M**2) + (0.166/M**4) + (0.035/M**6)
            C_D_base = 0.25/M
        C_D_stag = 0.85 * qRatio
        C_D_parasitic = max(1.3-0.3*(self.launchGuideLength/self.launchGuideOuterDiameter), 1) * C_D_stag
        outer_Area = np.pi * (self.launchGuideOuterDiameter/2)**2
        inner_Area = np.pi * (self.launchGuideInnerDiameter/2)**2
        A_parasitic = outer_Area - inner_Area
        if np.linalg.norm(airVelocity) == 0:
            return 0
        else:
            f = 0.25/((np.log10((self.launchGuideRoughness/(3.7*self.launchGuideInnerDiameter)) + (5.74/(self.reynoldsNumber(airVelocity, position, self.launchGuideInnerDiameter)**0.9))))**2)
        C_D_tube = f * (self.launchGuideLength/self.launchGuideInnerDiameter)
        return (C_D_tube * inner_Area + 0.7*(C_D_parasitic + C_D_base)*A_parasitic)/self.A_ref

    def axialDragCoefficient(self, time, airVelocity, position, AOA):
        if abs(np.degrees(AOA)) <= 17:
            scalefunction = (-3/24565)*abs(np.degrees(AOA))**3 + (9/2890)*(np.degrees(AOA))**2 + 1
        else:
            scalefunction = (13/1945085)*(np.degrees(AOA)-90)**3 + (2847/3890170)*(np.degrees(AOA)-90)**2
        return self.zeroAngleDragCoefficient(time, airVelocity, position) * scalefunction

    def axialDragForce(self, time, airVelocity, position, AOA):
        rho = self.airDensity(position)
        V = self.speed(airVelocity)
        if airVelocity[2] >= 0:
            A, C_A = self.A_ref, self.axialDragCoefficient(time, airVelocity, position, AOA)
        elif position[2] > self.MainDeploymentAltitude:
            A, C_A = self.parachuteArea, self.parachuteDragCoefficient
        else:
            A, C_A = self.MainArea, self.MainDragCoefficient
        return (1/2)*rho*V**2 * A * C_A

    # ------------------------------------------------------------------
    # Normal force / center of pressure
    # ------------------------------------------------------------------
    def normalForceFinCoefficientDerivative(self, airVelocity, position):
        M = self.machNumber(airVelocity, position)
        midchordSweepAngle = np.arctan((self.sweepLength + self.tipChord/2 - self.rootChord/2) / self.finSpan) if self.finSpan else 0.0
        beta = np.sqrt(abs(1-M**2))
        C_Nalpha0 = 2*np.pi/beta
        AR = 2*(self.finSpan**2)/self.A_fin
        F_D = AR/((1/(2*np.pi))*C_Nalpha0*np.cos(midchordSweepAngle))
        C_N_1fin = (C_Nalpha0 * F_D * (self.A_fin/self.A_ref) * np.cos(midchordSweepAngle))/(2+F_D*np.sqrt(1+4/(F_D**2)))    #(2*np.pi*((finSpan**2)/A_ref))/(2+((beta*finSpan**2)/(A_fin*np.cos(midchordSweepAngle))))
        C_N_fins = (self.numberOfFins/2) * C_N_1fin
        K_TB = 1 + (self.bodyRadius/(self.finSpan+self.bodyRadius))
        return K_TB * C_N_fins


    def normalForceCoefficientDerivative(self, airVelocity, position,AOA):
        C_N_finsWithInterference = self.normalForceFinCoefficientDerivative(airVelocity, position)
        C_N_nose = 2*np.cos(AOA)
        return C_N_nose + C_N_finsWithInterference

    def pitchMomentCoefficientDerivative(self, airVelocity, position, AOA):
        C_m_nose = (2*np.sin(AOA))/(self.A_ref*self.bodyDiameter)*(self.noseLength*self.A_ref-self.noseVolume)
        M = self.machNumber(airVelocity, position)
        AR = 2*self.finSpan**2/self.A_fin
        beta = np.sqrt(abs(1-M**2))

        A = np.array([
        [1,0.5,0.5**2,0.5**3,0.5**4,0.5**5],
        [0,1,2*0.5,3*0.5**2,4*0.5**3,5*0.5**4],
        [1,2,2**2,2**3,2**4,2**5],
        [0,1,2*2,3*2**2,4*2**3,5*2**4],
        [0,0,2,6*2,12*2**2, 20*2**3],
        [0,0,0,6,24*2,60*2**2]
        ])

        b = np.array([
            0.25, #p(0.5)
            0,  #p'(0.5)
            self.meanAerodynamicChordLengthOfFin*((AR*np.sqrt(abs(2**2-1)) - 0.67)/(2*AR*np.sqrt(abs(2**2-1))-1)),   #p(2) = f(2)
            ((0.68*AR)/(np.sqrt(3)*(2*AR*np.sqrt(3)-1)**2)),  #p'(2) = f'(2)
            0,  #p''(2)
            0   #p'''(2)
        ])

        coefficients = np.linalg.solve(A,b)

        rootChord, tipChord, sweepLength = self.rootChord, self.tipChord, self.sweepLength
        if M <= 0.5:
            x_f = (sweepLength/3)*((rootChord+2*tipChord)/(rootChord+tipChord))+(1/6)*((rootChord**2+tipChord**2+rootChord*tipChord)/(rootChord+tipChord))
        elif M <=2:
            x_f = coefficients[0] + coefficients[1]*M + coefficients[2]*M**2 + coefficients[3]*M**3 + coefficients[4]*M**4 + coefficients[5]*M**5
        else:
            x_f = self.meanAerodynamicChordLengthOfFin*((AR*beta - 0.67)/(2*AR*beta-1))
        x_fin = self.leadingEdgeToNosecone + x_f
        C_m_fins = self.normalForceFinCoefficientDerivative(airVelocity, position) * x_fin/self.bodyDiameter
        return C_m_nose+C_m_fins

    def normalForce(self, airVelocity, position, AOA):
        rho = self.airDensity(position)
        V = self.speed(airVelocity)
        C_N_alpha = self.normalForceCoefficientDerivative(airVelocity, position, AOA)
        return (1/2)*rho*V**2 * self.A_ref * C_N_alpha * AOA

    def centerOfPressure(self, airVelocity, position, AOA):
        if AOA != 0:
            return (self.pitchMomentCoefficientDerivative(airVelocity, position, AOA)/self.normalForceCoefficientDerivative(airVelocity, position, AOA))*self.bodyDiameter
        else:
            return self.CP

    # ------------------------------------------------------------------
    # Mass properties, thrust, launch rod
    # ------------------------------------------------------------------
    def moment_of_inertia(self, time):
        # Assuming a simple cylindrical rocket for now
        Mass = self.mass(time)
        CG = self.centerOfGravity(time)

        I_xx = Mass * CG * (self.rocketLength-CG)
        I_yy = I_xx
        I_zz = Mass * self.bodyRadius**2

        return ([I_xx, I_yy, I_zz])

    def on_rod(self, position):
        if self.leftRod:
            return False
        if np.linalg.norm(position) <= self.rodLength:
            return True
        self.leftRod = True
        return False

    def thrust_magnitude(self, time):
        return float(np.interp(time, self.timeCurve, self.thrustCurve, right=0.0))

    def mass(self, time):
        timeCurve, thrustCurve, cumulativeImpulse = self.timeCurve, self.thrustCurve, self.cumulativeImpulse
        if time <= timeCurve[0]:
            currentImpulse = self.thrust_magnitude(time) * time
        elif time >= timeCurve[-1]:
            currentImpulse = self.netImpulse
        else:
            i = np.searchsorted(timeCurve, time)
            if timeCurve[i] == time:
                currentImpulse = cumulativeImpulse[i]
            else:
                i -= 1
                currentImpulse = cumulativeImpulse[i] + 0.5 * (self.thrust_magnitude(time) + thrustCurve[i]) * (time - timeCurve[i])

        currentFuelMass = self.fuelMass * (1 - (currentImpulse / self.netImpulse))
        return self.dryMass + currentFuelMass

    def centerOfGravity(self, time):
        return (self.CG_dry*self.dryMass+self.CT*(self.mass(time)-self.dryMass))/(self.mass(time))

    # --- CP uncertainty (calibers) ---

    def cpSigma(self, M):
        return self.cpSigmaBase + self.cpSigmaPeak*np.exp(-((M-1)/self.cpWidth)**2)

    def cpOffset(self, time, M):
        """CP shift in metres (+ = aft = more stable)."""
        bump = np.exp(-((M-1)/self.cpWidth)**2)
        wander = float(np.interp(time, self.cpNoiseTimes, self.cpNoise)) if self.cpJitter else 0.0
        return (self.cpBiasZ*self.cpSigma(M) + self.cpJitter*bump*wander) * self.bodyDiameter

    # ------------------------------------------------------------------
    # Equations of motion
    # ------------------------------------------------------------------
    def derivatives(self, time, position, velocity, orientation, angular_velocity):
        wind = float(np.interp(time, self.timeWindSpeeds, self.windSpeeds, right=0.0))
        windVector = wind*np.array([np.sin(self.windHeading), np.cos(self.windHeading), 0])
        airFlow_velocity_world = velocity - windVector
        airFlow_velocity_body = inverse_rotate_vector(orientation, airFlow_velocity_world)
        AOA = self.angle_of_attack(airFlow_velocity_body)
        CP = self.centerOfPressure(airFlow_velocity_world, position, AOA) + self.cpOffset(time, self.machNumber(airFlow_velocity_world, position))
        CG = self.centerOfGravity(time)

        airFlow_speed = np.linalg.norm(airFlow_velocity_body)


        perpindicular_airFlow = airFlow_velocity_world - airFlow_speed * np.cos(AOA) * rotate_vector(orientation, [0,0,1])
        if AOA != 0:
            unit_perpindicular_airFlow = (1/np.linalg.norm(perpindicular_airFlow)) * perpindicular_airFlow
        else:
            unit_perpindicular_airFlow = np.zeros(3)

        normal_Force_body = inverse_rotate_vector(orientation, -self.normalForce(airFlow_velocity_world, position, AOA)*unit_perpindicular_airFlow)

        axialDrag_Force = self.axialDragForce(time, airFlow_velocity_world, position, AOA)
        if airFlow_velocity_world[2] < 0 and airFlow_speed > 0:
            aero_force_world = -axialDrag_Force * airFlow_velocity_world / airFlow_speed
            aero_force_body = inverse_rotate_vector(orientation, aero_force_world)
            moment_scale = 0.0
        else:
            aero_force_body = normal_Force_body + (np.array([0,0,axialDrag_Force]) * -np.sign(airFlow_velocity_body[2]))
            aero_force_world = rotate_vector(orientation, aero_force_body)
            moment_scale = 1.0
        thrust_body = np.array([0.0, 0.0, self.thrust_magnitude(time)])
        thrust_world = rotate_vector(orientation, thrust_body)

        gravity_world = np.array([0.0, 0.0, -self.mass(time) * self.gravity])
        if self.on_rod(position):
            netForce_world = thrust_world
        else:
            netForce_world = thrust_world + gravity_world + aero_force_world

        # Acceleration
        acceleration = netForce_world / self.mass(time)

        # Moments
        total_moment_body = moment_scale * np.cross([0,0,CG-CP], aero_force_body)

        # Angular acceleration
        I_matrix = np.diag(self.moment_of_inertia(time))
        angular_acceleration = np.linalg.solve(I_matrix, total_moment_body - np.cross(angular_velocity, I_matrix @ angular_velocity))

        # Quaternion derivative
        omega_quaternion = np.array([
            0.0,
            angular_velocity[0],
            angular_velocity[1],
            angular_velocity[2]
        ])

        ori_dt = 0.5 * quaternion_multiply(orientation,omega_quaternion)

        return (velocity, acceleration, ori_dt, angular_acceleration)

    def rk4_step(self, time, position, velocity, orientation, angular_velocity):
        timeStep = self.timeStep
        derivatives = self.derivatives

        k1_pos, k1_vel, k1_q, k1_omega = derivatives(time, position, velocity, orientation, angular_velocity)

        k2_pos, k2_vel, k2_q, k2_omega = derivatives(time + timeStep/2, position + k1_pos * timeStep/2, velocity + k1_vel * timeStep/2, orientation + k1_q * timeStep/2, angular_velocity + k1_omega * timeStep/2)

        k3_pos, k3_vel, k3_q, k3_omega = derivatives(time + timeStep/2, position + k2_pos * timeStep/2, velocity + k2_vel * timeStep/2, orientation + k2_q * timeStep/2, angular_velocity + k2_omega * timeStep/2)

        k4_pos, k4_vel, k4_q, k4_omega = derivatives(time + timeStep, position + k3_pos * timeStep, velocity + k3_vel * timeStep, orientation + k3_q * timeStep, angular_velocity + k3_omega * timeStep)


        new_position = position + (timeStep/6) * (k1_pos + 2*k2_pos + 2*k3_pos + k4_pos)
        new_velocity = velocity + (timeStep/6) * (k1_vel + 2*k2_vel + 2*k3_vel + k4_vel)
        new_orientation = orientation + (timeStep/6) * (k1_q + 2*k2_q + 2*k3_q + k4_q)
        new_angular_velocity = angular_velocity + (timeStep/6) * (k1_omega + 2*k2_omega + 2*k3_omega + k4_omega)

        new_orientation = quaternion_normalize(new_orientation)

        if(self.on_rod(new_position)):
                new_orientation = orientation
                new_angular_velocity = np.zeros(3)

        return (new_position,new_velocity,new_orientation,new_angular_velocity)

    def run_simulation(self):
        # The web app runs many times per instance, so start every run from
        # the launch pad with the current parameters and fresh wind.
        self.calculate_values()
        self.generate_wind()
        self.generate_cp_noise()
        self.reset_state()

        time, position, velocity = self.time, self.position, self.velocity
        orientation, angular_velocity = self.orientation, self.angular_velocity
        tempNum = 0
        while -0.01 <= position[2] <= 10000000:  #limits the sim to actual atmosphere, avoid errors.

            euler_orientation = quaternion_to_euler(orientation)
            wind = float(np.interp(time, self.timeWindSpeeds, self.windSpeeds, right=0.0))
            windVector = wind* np.array([np.sin(self.windHeading),np.cos(self.windHeading), 0])
            derivative = self.derivatives(time, position, velocity, orientation, angular_velocity)
            self.timeList.append(time)
            self.massList.append(self.mass(time))
            self.positionHistory = np.hstack((self.positionHistory, np.array(position).reshape(-1, 1)))
            self.velocityHistory = np.hstack((self.velocityHistory, np.array(velocity).reshape(-1, 1)))
            self.forceHistory = np.hstack((self.forceHistory, np.array(derivative[1] * self.mass(time)).reshape(-1, 1)))
            self.orientationHistory = np.hstack((self.orientationHistory, np.array(euler_orientation).reshape(-1,1)))
            self.angularVelocityHistory = np.hstack((self.angularVelocityHistory, np.array(angular_velocity).reshape(-1, 1)))
            self.momentHistory = np.hstack((self.momentHistory, np.array(derivative[3] * np.array(self.moment_of_inertia(time))).reshape(-1, 1)))
            self.machList.append(self.machNumber(velocity - windVector, position))
            self.AOAList.append(self.angle_of_attack(velocity))
            if velocity[2] >= 0:
                self.CDList.append(self.axialDragCoefficient(time, velocity - windVector, position, self.angle_of_attack(velocity)))
            else:
                self.CDList.append(self.parachuteDragCoefficient)

            position, velocity, orientation, angular_velocity = self.rk4_step(time, position, velocity, orientation, angular_velocity)

            time += self.timeStep
        for i in range(len(self.positionHistory[0])):
            self.planarHistory.append(np.sqrt(self.positionHistory[0][i]**2 + self.positionHistory[1][i]**2))
        self.time, self.position, self.velocity = time, position, velocity
        self.orientation, self.angular_velocity = orientation, angular_velocity

    # ------------------------------------------------------------------
    # Results for the frontend - read straight from the stored histories
    # ------------------------------------------------------------------
    # Column name: (unit, how it is read from the stored histories)
    RESULT_UNITS = {
        "time": "s",
        "mass": "kg",
        "pos_north": "m", "pos_east": "m", "pos_up": "m",
        "vel_north": "m/s", "vel_east": "m/s", "vel_up": "m/s",
        "speed": "m/s",
        "force_north": "N", "force_east": "N", "force_up": "N",
        "pitch": "deg", "yaw": "deg", "roll": "deg",
        "pitch_rate": "rad/s", "yaw_rate": "rad/s", "roll_rate": "rad/s",
        "pitch_moment": "N*m", "yaw_moment": "N*m", "roll_moment": "N*m",
        "mach": "-",
        "aoa": "deg",
        "cd": "-",
    }

    # Old 2D names main.py still uses (y_pos for the altitude chart)
    LEGACY_ALIASES = {"x_pos": "pos_north", "y_pos": "pos_up"}

    def get_results(self, include_legacy=True):
        """All stored time series as a dict of lists, one entry per time
        step. Units are in RESULT_UNITS. World frame is North, East, Up;
        rates and moments are body frame, as in src/sim.py."""
        p, v, F = self.positionHistory, self.velocityHistory, self.forceHistory
        o, w, M = self.orientationHistory, self.angularVelocityHistory, self.momentHistory
        results = {
            "time": list(self.timeList),
            "mass": list(self.massList),
            "pos_north": p[0].tolist(), "pos_east": p[1].tolist(), "pos_up": p[2].tolist(), "planar_distance": self.planarHistory,
            "vel_north": v[0].tolist(), "vel_east": v[1].tolist(), "vel_up": v[2].tolist(),
            "speed": np.linalg.norm(v, axis=0).tolist(),
            "force_north": F[0].tolist(), "force_east": F[1].tolist(), "force_up": F[2].tolist(),
            "pitch": o[0].tolist(), "yaw": o[1].tolist(), "roll": o[2].tolist(),
            "pitch_rate": w[0].tolist(), "yaw_rate": w[1].tolist(), "roll_rate": w[2].tolist(),
            "pitch_moment": M[0].tolist(), "yaw_moment": M[1].tolist(), "roll_moment": M[2].tolist(),
            "mach": [float(x) for x in self.machList],
            "aoa": np.degrees(self.AOAList).tolist(),
            "cd": [float(x) for x in self.CDList],
        }
        if include_legacy:
            for old, new in self.LEGACY_ALIASES.items():
                results[old] = list(results[new])
        return results

    def results_csv_bytes(self):
        """CSV of the last run's 3D data, built in memory."""
        results = self.get_results(include_legacy=False)
        buf = io.StringIO()
        writer = csv.writer(buf)
        writer.writerow(results.keys())
        writer.writerows(zip(*results.values()))
        return buf.getvalue().encode("utf-8")

    def get_summary(self):
        """The same figures src/sim.py prints at the end of a run."""
        if not self.timeList:
            raise RuntimeError("run_simulation() must be called before get_summary()")
        p, v = self.positionHistory, self.velocityHistory

        rodVelocity = 0
        for i in range(len(p[0])-1):
            dist1=np.linalg.norm(p[:, i])
            dist2=np.linalg.norm(p[:, i+1])
            if dist1 < self.rodLength <= dist2:
                rodVelocity = np.linalg.norm((v[:, i]+v[:, i+1])/2)
                break

        apogee_index = int(np.argmax(p[2]))
        return {
            "rod_velocity": float(rodVelocity),
            "apogee": float(max(p[2])),
            "max_velocity": float(max(v[2])),
            "max_acceleration": float(max(self.forceHistory[2] / np.array(self.massList))),
            "time_to_apogee": float(self.timeList[apogee_index]),
            "flight_time": float(self.timeList[-1]),
            "ground_hit_velocity": float(v[2][-1]),
            "range": float(np.sqrt(p[0][-1]**2 + p[1][-1]**2)),
            "max_mach": float(max(self.machList)),
        }