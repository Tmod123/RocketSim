#sim.py - V3.0 (class-based, for per-session isolation)
#
# Previously this module stored all simulation state (dryMass, timeCurve,
# result lists, etc.) as module-level globals. In Streamlit, a module is
# imported once per server process and shared by every visitor, so one
# user's parameter changes / simulation runs were overwriting everyone
# else's. Wrapping the state in a class lets main.py create one instance
# per browser session (via st.session_state) so sessions never touch each
# other's data. It also means we never need to write motor/result files
# to a shared path on disk - motor data and results just live in memory.

import numpy as np
import csv
import io


class RocketSim:

    def __init__(self, motor_source="Motors/Estes_C6.csv", **overrides):
        # ---- Constants ----
        self.gravity = 9.81  # m/s^2
        self.timeStep = 0.05  # seconds

        # ---- Rocket geometry (not exposed to the GUI) ----
        self.L = 0.10
        self.NBD = 0.025
        self.d11 = 0.0
        self.d21 = 0.0
        self.d12 = 0.0
        self.d22 = 0.0
        self.Lt1 = 0.0
        self.Lt2 = 0.0
        self.xs1 = 0.0
        self.xs2 = 0.0
        self.a1 = 0.0508
        self.b1 = 0.0508
        self.a2 = 0.0
        self.b2 = 0.0
        self.f1 = 3.0
        self.f2 = 0.0
        self.Y1 = 0.03931
        self.Y2 = 0.0
        self.R1 = 0.0125
        self.R2 = 0.0
        self.Xf1 = 0.3492
        self.Xf2 = 0.0
        self.m1 = 0.0254
        self.m2 = 0.0
        self.s1 = 0.03
        self.s2 = 0.0
        self.finSets = 1
        self.transitions = 0
        self.noseType = "ogive"  # cone, ogive, paraboloid, ellipsoid

        # ---- Editable rocket parameters ----
        self.dryMass = 0.0605
        self.fuelMass = 0.011
        self.radius = 0.0127
        self.dragCoefficient = 0.634
        self.parachuteArea = 0.0707
        self.parachuteDragCoefficient = 0.80
        self.rodLength = 1.0
        self.launchAngle = np.radians(5)
        self.windSpeed = 0
        self.rocketLength = 0.425
        self.CG_dry = 0.24
        self.CG_wet = 0.26
        self.CN_alpha = 11.97

        for key, value in overrides.items():
            setattr(self, key, value)

        self.CP = self.center_of_pressure()

        self.timeCurve, self.thrustCurve = [], []
        self.set_motor_data(motor_source)

        self._reset_state()

    # ------------------------------------------------------------------
    # Motor data
    # ------------------------------------------------------------------
    @staticmethod
    def _parse_motor_data(file_obj):
        times = [0.000]
        thrusts = [0.000]
        for row in csv.reader(file_obj):
            if not row:
                continue
            try:
                times.append(float(row[0]))
                thrusts.append(float(row[1]))
            except (ValueError, IndexError):
                pass
        return times, thrusts

    def set_motor_data(self, source):
        """Load a thrust curve. `source` can be a filesystem path (str) or an
        in-memory text stream (e.g. io.StringIO built from an uploaded file's
        bytes) - so an uploaded motor file never has to touch disk."""
        if isinstance(source, str):
            with open(source, newline="") as f:
                self.timeCurve, self.thrustCurve = self._parse_motor_data(f)
        else:
            source.seek(0)
            self.timeCurve, self.thrustCurve = self._parse_motor_data(source)

    # ------------------------------------------------------------------
    # Aerodynamics setup
    # ------------------------------------------------------------------
    def center_of_pressure(self):
        CN_n = 2
        if self.finSets >= 1:
            if self.finSets == 2:
                CN_f2 = (1 + self.R2 / (self.s2 + self.R2)) * (
                    (4 * self.f2 * (self.s2 / self.NBD)**2) /
                    (1 + np.sqrt(1 + ((2 * self.Y2) /
                                      (self.a2 + self.b2))**2)))
                Pf2 = CN_f2 * (self.Xf2 + (self.m2 * (self.a2 + 2 * self.b2) /
                                           (3 *
                                            (self.a2 + self.b2))) + (1 / 6) *
                               (self.a2 + self.b2 - (self.a2 * self.b2) /
                                (self.a2 + self.b2)))
            else:
                CN_f2, Pf2 = 0, 0
            CN_f1 = (1 + self.R1 / (self.s1 + self.R1)) * (
                (4 * self.f1 * (self.s1 / self.NBD)**2) /
                (1 + np.sqrt(1 + ((2 * self.Y1) / (self.a1 + self.b1))**2)))
            Pf1 = CN_f1 * (self.Xf1 + (self.m1 * (self.a1 + 2 * self.b1) /
                                       (3 * (self.a1 + self.b1))) + (1 / 6) *
                           (self.a1 + self.b1 - (self.a1 * self.b1) /
                            (self.a1 + self.b1)))
        else:
            CN_f1, Pf1, CN_f2, Pf2 = 0, 0, 0, 0

        if self.transitions >= 1:
            if self.transitions == 2:
                CN_s2 = 2 * ((self.d22 / self.NBD)**2 -
                             (self.d12 / self.NBD)**2)
                Ps2 = CN_s2 * (self.xs2 + (self.Lt2 / 3) *
                               (1 + (1 - (self.d12 / self.d22)) /
                                (1 - (self.d12 / self.d22)**2)))
            else:
                CN_s2, Ps2 = 0, 0
            CN_s1 = 2 * ((self.d21 / self.NBD)**2 - (self.d11 / self.NBD)**2)
            Ps1 = CN_s1 * (self.xs1 + (self.Lt1 / 3) *
                           (1 + (1 - (self.d11 / self.d21)) /
                            (1 - (self.d11 / self.d21)**2)))
        else:
            CN_s1, Ps1, CN_s2, Ps2 = 0, 0, 0, 0

        if self.noseType == "cone":
            Pn = CN_n * (0.6667 * self.L)
        elif self.noseType == "ogive":
            Pn = CN_n * (0.466 * self.L)
        elif self.noseType == "paraboloid":
            Pn = CN_n * (0.5 * self.L)
        elif self.noseType == "ellipsoid":
            Pn = CN_n * (0.3333 * self.L)
        else:
            Pn = 0

        netNormalForce = CN_n + CN_s1 + CN_s2 + CN_f1 + CN_f2
        netMoment = Pn + Ps1 + Ps2 + Pf1 + Pf2
        return netMoment / netNormalForce

    # ------------------------------------------------------------------
    # Integrator state
    # ------------------------------------------------------------------
    def _reset_state(self):
        self.time = 0
        self.y_pos = 0
        self.y_velocity = 0
        self.x_pos = 0
        self.x_velocity = 0
        self.theta = self.launchAngle
        self.omega = 0
        self.leftRod = False

        self.timeList = []
        self.massList = []
        self.x_posList = []
        self.y_posList = []
        self.x_velocityList = []
        self.y_velocityList = []
        self.thetaList = []
        self.omegaList = []
        self.alphaList = []
        self.torqueList = []
        self.inertiaList = []
        self.normalForceList = []
        self.x_accelerationList = []
        self.y_accelerationList = []
        self.x_netForceList = []
        self.y_netForceList = []
        self.weightList = []
        self.x_thrustList = []
        self.y_thrustList = []
        self.x_dragList = []
        self.y_dragList = []
        self.densityList = []

    # ------------------------------------------------------------------
    # Physics
    # ------------------------------------------------------------------
    def moment_of_inertia(self, current_mass, time):
        I_center = (1 / 12) * current_mass * self.rocketLength**2
        d = self.center_of_gravity(time) - self.rocketLength / 2
        return I_center + current_mass * d**2

    def angle_of_attack(self, theta, x_velocity, y_velocity):
        return np.arctan2(x_velocity, y_velocity) - theta

    def normal_force(self, theta, y_pos, x_velocity, y_velocity):
        dynamicPressure = 0.5 * self.air_density(y_pos) * self.net_speed(
            x_velocity, y_velocity)**2
        Area = np.pi * self.radius**2
        return dynamicPressure * Area * self.CN_alpha * self.angle_of_attack(
            theta, x_velocity, y_velocity)

    def yaw_torque(self, theta, y_pos, x_velocity, y_velocity, time):
        return self.normal_force(theta, y_pos, x_velocity, y_velocity) * (
            self.CP - self.center_of_gravity(time))

    def alpha(self, time, theta, y_pos, x_velocity, y_velocity):
        return self.yaw_torque(theta, y_pos, x_velocity, y_velocity,
                               time) / self.moment_of_inertia(
                                   self.mass(time), time)

    def air_density(self, y_pos):
        return 1.225 * (1 -
                        ((0.0065 / 288.15) * y_pos))**((self.gravity /
                                                        (287.05 * 0.0065)) - 1)

    def net_speed(self, x_velocity, y_velocity):
        return np.sqrt((x_velocity**2) + (y_velocity**2))

    def on_rod(self, x_pos, y_pos):
        if self.leftRod:
            return False
        if np.sqrt((x_pos**2) + (y_pos**2)) <= self.rodLength:
            return True
        self.leftRod = True
        return False

    def thrust_magnitude(self, time):
        Output = 0
        
        if time > self.timeCurve[-1]:
            return 0
        if time <= self.timeCurve[0]:
            print(self.thrustCurve[0])
            print(self.timeCurve[0])
            try:
                Output = self.thrustCurve[1] * time / self.timeCurve[1]
                
            except ZeroDivisionError:
                Output = 0
                
            return Output
        for i in range(len(self.timeCurve) - 1):
            if self.timeCurve[i] < time <= self.timeCurve[i + 1]:
                Time1, Time2 = self.timeCurve[i], self.timeCurve[i + 1]
                Thrust1, Thrust2 = self.thrustCurve[i], self.thrustCurve[i + 1]
                return Thrust1 + (Thrust2 -
                                  Thrust1) * (time - Time1) / (Time2 - Time1)
        return 0

    def x_component(self, magnitude, theta):
        return magnitude * np.sin(theta)

    def y_component(self, magnitude, theta):
        return magnitude * np.cos(theta)

    def acceleration_x(self, time, y_pos, x_velocity, y_velocity, theta):
        return self.net_force_x(time, y_pos, x_velocity, y_velocity,
                                theta) / self.mass(time)

    def acceleration_y(self, time, y_pos, x_velocity, y_velocity, theta):
        return self.net_force_y(time, y_pos, x_velocity, y_velocity,
                                theta) / self.mass(time)

    def net_force_x(self, time, y_pos, x_velocity, y_velocity, theta):
        if y_pos == 0:
            return self.x_component(self.thrust_magnitude(time), theta)
        return self.x_component(
            self.thrust_magnitude(time), theta) - self.x_component(
                self.drag_magnitude(y_pos, x_velocity, y_velocity),
                np.arctan2(x_velocity, y_velocity))

    def net_force_y(self, time, y_pos, x_velocity, y_velocity, theta):
        if y_pos == 0:
            return self.y_component(self.thrust_magnitude(time), theta)
        return (self.y_component(self.thrust_magnitude(time), theta) -
                self.weight(time) - self.y_component(
                    self.drag_magnitude(y_pos, x_velocity, y_velocity),
                    np.arctan2(x_velocity, y_velocity)))

    def drag_magnitude(self, y_pos, x_velocity, y_velocity):
        speed = self.net_speed(x_velocity, y_velocity)
        if y_velocity >= 0:
            A, cd = 3.14 * self.radius**2, self.dragCoefficient
        else:
            A, cd = self.parachuteArea, self.parachuteDragCoefficient
        return 0.5 * self.air_density(y_pos) * cd * A * speed**2

    def weight(self, time):
        return self.mass(time) * self.gravity

    def mass(self, time):
        currentImpulse = 0
        netImpulse = 0
        for i in range(1, len(self.timeCurve)):
            if time < self.timeCurve[0]:
                currentImpulse = self.thrust_magnitude(time) * time
            elif self.timeCurve[i] <= time:
                currentImpulse += 0.5 * (
                    self.thrustCurve[i] + self.thrustCurve[i - 1]) * (
                        self.timeCurve[i] - self.timeCurve[i - 1])
            else:
                currentImpulse += 0.5 * (self.thrust_magnitude(time) +
                                         self.thrustCurve[i - 1]) * (
                                             time - self.timeCurve[i - 1])
                break
        for i in range(1, len(self.timeCurve)):
            netImpulse += 0.5 * (self.thrustCurve[i] +
                                 self.thrustCurve[i - 1]) * (
                                     self.timeCurve[i] - self.timeCurve[i - 1])
        currentFuelMass = self.fuelMass * (1 - (currentImpulse / netImpulse))
        return self.dryMass + currentFuelMass

    def center_of_gravity(self, time):
        # CT is recomputed here (not cached) so edits to dryMass/fuelMass/CG_dry/CG_wet
        # from the sidebar take effect on the next run_simulation() call.
        CT = (self.CG_wet * (self.dryMass + self.fuelMass) -
              self.CG_dry * self.dryMass) / self.fuelMass
        m = self.mass(time)
        return (self.CG_dry * self.dryMass + CT * (m - self.dryMass)) / m

    # ------------------------------------------------------------------
    # Integration
    # ------------------------------------------------------------------
    def rk4_step(self, time, x_pos, y_pos, x_velocity, y_velocity, theta,
                 omega):
        ts = self.timeStep

        k1_x = x_velocity
        k1_y = y_velocity
        k1_vx = self.acceleration_x(time, y_pos, x_velocity, y_velocity, theta)
        k1_vy = self.acceleration_y(time, y_pos, x_velocity, y_velocity, theta)
        k1_0 = omega
        k1_w = self.alpha(time, theta, y_pos, x_velocity, y_velocity)

        k2_x = x_velocity + k1_vx * ts / 2
        k2_y = y_velocity + k1_vy * ts / 2
        k2_vx = self.acceleration_x(time + ts / 2, y_pos + k1_y * ts / 2,
                                    x_velocity + k1_vx * ts / 2,
                                    y_velocity + k1_vy * ts / 2,
                                    theta + k1_0 * ts / 2)
        k2_vy = self.acceleration_y(time + ts / 2, y_pos + k1_y * ts / 2,
                                    x_velocity + k1_vx * ts / 2,
                                    y_velocity + k1_vy * ts / 2,
                                    theta + k1_0 * ts / 2)
        k2_0 = omega + k1_w * ts / 2
        k2_w = self.alpha(time + ts / 2, theta + k1_0 * ts / 2,
                          y_pos + k1_y * ts / 2, x_velocity + k1_vx * ts / 2,
                          y_velocity + k1_vy * ts / 2)

        k3_x = x_velocity + k2_vx * ts / 2
        k3_y = y_velocity + k2_vy * ts / 2
        k3_vx = self.acceleration_x(time + ts / 2, y_pos + k2_y * ts / 2,
                                    x_velocity + k2_vx * ts / 2,
                                    y_velocity + k2_vy * ts / 2,
                                    theta + k2_0 * ts / 2)
        k3_vy = self.acceleration_y(time + ts / 2, y_pos + k2_y * ts / 2,
                                    x_velocity + k2_vx * ts / 2,
                                    y_velocity + k2_vy * ts / 2,
                                    theta + k2_0 * ts / 2)
        k3_0 = omega + k2_w * ts / 2
        k3_w = self.alpha(time + ts / 2, theta + k2_0 * ts / 2,
                          y_pos + k2_y * ts / 2, x_velocity + k2_vx * ts / 2,
                          y_velocity + k2_vy * ts / 2)

        k4_x = x_velocity + k3_vx * ts
        k4_y = y_velocity + k3_vy * ts
        k4_vx = self.acceleration_x(time + ts, y_pos + k3_y * ts,
                                    x_velocity + k3_vx * ts,
                                    y_velocity + k3_vy * ts, theta + k3_0 * ts)
        k4_vy = self.acceleration_y(time + ts, y_pos + k3_y * ts,
                                    x_velocity + k3_vx * ts,
                                    y_velocity + k3_vy * ts, theta + k3_0 * ts)
        k4_0 = omega + k3_w * ts
        k4_w = self.alpha(time + ts, theta + k3_0 * ts, y_pos + k3_y * ts,
                          x_velocity + k3_vx * ts, y_velocity + k3_vy * ts)

        new_x_pos = x_pos + (ts / 6) * (k1_x + 2 * k2_x + 2 * k3_x + k4_x)
        new_y_pos = y_pos + (ts / 6) * (k1_y + 2 * k2_y + 2 * k3_y + k4_y)
        new_x_velocity = x_velocity + (ts / 6) * (k1_vx + 2 * k2_vx +
                                                  2 * k3_vx + k4_vx)
        new_y_velocity = y_velocity + (ts / 6) * (k1_vy + 2 * k2_vy +
                                                  2 * k3_vy + k4_vy)
        new_theta = theta + (ts / 6) * (k1_0 + 2 * k2_0 + 2 * k3_0 + k4_0)
        new_omega = omega + (ts / 6) * (k1_w + 2 * k2_w + 2 * k3_w + k4_w)

        if self.on_rod(new_x_pos, new_y_pos):
            new_theta = theta
            new_omega = 0

        if new_y_velocity < 0:
            new_theta = theta
            new_omega = 0

        new_theta = (new_theta + np.pi) % (2 * np.pi) - np.pi

        return new_x_pos, new_y_pos, new_x_velocity, new_y_velocity, new_theta, new_omega

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------
    def run_simulation(self):
        self._reset_state()

        time, x_pos, y_pos = self.time, self.x_pos, self.y_pos
        x_velocity, y_velocity = self.x_velocity, self.y_velocity
        theta, omega = self.theta, self.omega

        while y_pos >= -0.01:
            current_mass = self.mass(time)
            current_thrust = self.thrust_magnitude(time)
            current_thrust_x = self.x_component(current_thrust, theta)
            current_thrust_y = self.y_component(current_thrust, theta)
            current_weight = current_mass * self.gravity
            current_drag = self.drag_magnitude(y_pos, x_velocity, y_velocity)
            current_drag_x = self.x_component(
                current_drag, np.arctan2(x_velocity, y_velocity))
            current_drag_y = self.y_component(
                current_drag, np.arctan2(x_velocity, y_velocity))
            current_normalForce = self.normal_force(theta, y_pos, x_velocity,
                                                    y_velocity)
            current_inertia = self.moment_of_inertia(
                current_mass, self.center_of_gravity(time))

            if y_pos == 0:
                current_netForce_x = current_thrust_x
                current_netForce_y = current_thrust_y
            else:
                current_netForce_x = current_thrust_x - current_drag_x
                current_netForce_y = current_thrust_y - current_weight - current_drag_y

            if self.on_rod(x_pos, y_pos):
                current_torque = 0
            else:
                current_torque = self.yaw_torque(theta, y_pos, x_velocity,
                                                 y_velocity, time)

            current_acceleration_x = current_netForce_x / current_mass
            current_acceleration_y = current_netForce_y / current_mass
            current_alpha = current_torque / current_inertia

            self.timeList.append(time)
            self.x_posList.append(x_pos)
            self.y_posList.append(y_pos)
            self.x_velocityList.append(x_velocity)
            self.y_velocityList.append(y_velocity)
            self.thetaList.append(np.degrees(theta))
            self.omegaList.append(omega)
            self.alphaList.append(current_alpha)
            self.torqueList.append(current_torque)
            self.inertiaList.append(current_inertia)
            self.normalForceList.append(current_normalForce)
            self.massList.append(current_mass)
            self.x_accelerationList.append(current_acceleration_x)
            self.y_accelerationList.append(current_acceleration_y)
            self.x_netForceList.append(current_netForce_x)
            self.y_netForceList.append(current_netForce_y)
            self.weightList.append(current_weight)
            self.x_thrustList.append(current_thrust_x)
            self.y_thrustList.append(current_thrust_y)
            self.x_dragList.append(current_drag_x)
            self.y_dragList.append(current_drag_y)
            self.densityList.append(self.air_density(y_pos))

            x_pos, y_pos, x_velocity, y_velocity, theta, omega = self.rk4_step(
                time, x_pos, y_pos, x_velocity, y_velocity, theta, omega)
            time += self.timeStep

        self.time, self.x_pos, self.y_pos = time, x_pos, y_pos
        self.x_velocity, self.y_velocity = x_velocity, y_velocity
        self.theta, self.omega = theta, omega

    def get_results(self):
        """Return all recorded time-series data as a dict of lists (used to
        build a pandas DataFrame in main.py - no shared file needed)."""
        return {
            "time": self.timeList,
            "mass": self.massList,
            "x_pos": self.x_posList,
            "y_pos": self.y_posList,
            "x_vel": self.x_velocityList,
            "velocity": self.y_velocityList,
            "x_acc": self.x_accelerationList,
            "acceleration": self.y_accelerationList,
            "x_net": self.x_netForceList,
            "net force": self.y_netForceList,
            "weight": self.weightList,
            "x_thrust": self.x_thrustList,
            "thrust": self.y_thrustList,
            "x_drag": self.x_dragList,
            "drag": self.y_dragList,
            "theta": self.thetaList,
            "omega": self.omegaList,
            "alpha": self.alphaList,
            "torque": self.torqueList,
            "inertia": self.inertiaList,
            "normalF": self.normalForceList,
            "air density": self.densityList,
        }

    def results_csv_bytes(self):
        """Build a CSV of the last run's results entirely in memory (no shared
        output.csv on disk, so one session's download never races another's)."""
        results = self.get_results()
        buf = io.StringIO()
        writer = csv.writer(buf)
        writer.writerow(results.keys())
        writer.writerows(zip(*results.values()))
        return buf.getvalue().encode("utf-8")

    def get_summary(self):
        if not self.timeList:
            raise RuntimeError(
                "run_simulation() must be called before get_summary()")

        rodVelocity = 0
        for i in range(len(self.y_posList) - 1):
            dist1 = np.sqrt(self.y_posList[i]**2 + self.x_posList[i]**2)
            dist2 = np.sqrt(self.y_posList[i + 1]**2 +
                            self.x_posList[i + 1]**2)
            if dist1 < self.rodLength <= dist2:
                rodVelocity = np.sqrt(
                    (self.y_velocityList[i]**2 +
                     2 * self.y_accelerationList[i] *
                     (self.rodLength * np.cos(self.launchAngle) -
                      self.y_posList[i])) +
                    (self.x_velocityList[i]**2 +
                     2 * self.x_accelerationList[i] *
                     (self.rodLength * np.sin(self.launchAngle) -
                      self.x_posList[i])))
                break

        apogee = max(self.y_posList)
        return {
            "rod_velocity": rodVelocity,
            "apogee": apogee,
            "max_velocity": max(self.y_velocityList),
            "max_acceleration": max(self.y_accelerationList),
            "time_to_apogee": self.timeList[self.y_posList.index(apogee)],
            "flight_time": self.timeList[-1],
            "ground_hit_velocity": self.y_velocityList[-1],
            "range": self.x_posList[-1],
        }
