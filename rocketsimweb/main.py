import io #agar
import plots
import streamlit
import pandas
import numpy as np
import matplotlib.pyplot as pyplot
from threading import RLock

import sim

streamlit.set_page_config(page_title="rocket flight simulator", layout="wide")
streamlit.title("rocket flight simulator")

# One RocketSim instance per browser session. st.session_state is scoped to a
# single user's session (unlike a plain module-level variable, which the
# whole server process - i.e. every visitor - would share).
if "sim" not in streamlit.session_state:
    streamlit.session_state.sim = sim.RocketSim()
    


def run_button_pressed():
    rocket = streamlit.session_state.sim
    rocket.run_simulation()
    summary = rocket.get_summary()
    streamlit.session_state.summary = summary

    # Keep results in memory (session_state), not a shared file on disk.
    streamlit.session_state.results = rocket.get_results()


def update_variable(variable_name, state_key):
    setattr(streamlit.session_state.sim, variable_name, streamlit.session_state[state_key])

summary_container = streamlit.container(border=True)
if "summary" in streamlit.session_state:
    s = streamlit.session_state.summary
    summary_container.markdown(
        f"<u>**Rod velocity**</u> is {round(s['rod_velocity'], 3)}m/s<br>"
        f"<u>**Apogee**</u> is at {round(s['apogee'], 3)} meters<br>"
        f"<u>**Time to apogee**</u> is {round(s['time_to_apogee'], 3)} seconds<br>"
        f"<u>**Max velocity**</u> is {round(s['max_velocity'], 3)}m/s<br>"
        f"<u>**Max acceleration**</u> is {round(s['max_acceleration'], 3)}m/s<sup>2</sup><br>"
        f"<u>**Flight time**</u> is {round(s['flight_time'], 3)} seconds<br>"
        f"<u>**Final velocity**</u> is {round(s['ground_hit_velocity'], 3)} m/s<br>"
        f"**Max mach was <u>{s['max_mach']}</u>**<br>"
        f"**Rocket landed {round(s['range'])} meters from the pad**",
        
        True)


with streamlit.sidebar:
    streamlit.button("run simulation", on_click=run_button_pressed)

    uploaded_file = streamlit.file_uploader("Choose a motor file", type="csv")
    if uploaded_file is not None:
        # Parse straight from the uploaded bytes - never touches disk, so
        # concurrent users can't clobber each other's motor file.
        text = uploaded_file.getvalue().decode("utf-8")
        streamlit.session_state.sim.set_motor_data(io.StringIO(text))

    dry_mass = streamlit.number_input(
        "parachute drag coefficient",
        value=0.8,
        format="%.5f",
        key="parachute_cd_input",
        on_change=update_variable,
        args=("parachuteDragCoefficient", "parachute_cd_input")
    )
    dry_mass = streamlit.number_input(
        "parachute area (m^2)",
        value=0.07306,
        format="%.5f",
        key="parachute_area_input",
        on_change=update_variable,
        args=("parachuteArea", "parachute_area_input")
    )
    dry_mass = streamlit.number_input(
        "main deployment altitude (m)",
        value=305,
        format="%.5f",
        key="main_deployment_altitude_input",
        on_change=update_variable,
        args=("MainDeploymentAltitude", "main_deployment_altitude_input")
    )
    dry_mass = streamlit.number_input(
        "main drag coefficient",
        value=1.550,
        format="%.5f",
        key="main_drag_coefficient_input",
        on_change=update_variable,
        args=("MainDragCoefficient", "main_drag_coefficient_input")
    )
    dry_mass = streamlit.number_input(
        "main area (m^2)",
        value=1.169,
        format="%.5f",
        key="main_area_input",
        on_change=update_variable,
        args=("MainArea", "main_area_input")
    )
    dry_mass = streamlit.number_input(
        "surface roughness (m)",
        value=0.00002,
        format="%.5f",
        key="surface_roughness_input",
        on_change=update_variable,
        args=("surfaceRoughness", "surface_roughness_input")
    )
    dry_mass = streamlit.number_input(
        "launch guide roughness",
        value=0.00006,
        format="%.5f",
        key="launch_guide_roughness_input",
        on_change=update_variable,
        args=("launchGuideRoughness", "launch_guide_roughness_input")
    )
    dry_mass = streamlit.number_input(
        "dry mass (kg)",
        value=4.991,
        format="%.5f",
        key="dry_mass_input",
        on_change=update_variable,
        args=("dryMass", "dry_mass_input")
    )
    dry_mass = streamlit.number_input(
        "fuel_mass (kg)",
        value=1.552,
        format="%.5f",
        key="fuel_mass_input",
        on_change=update_variable,
        args=("fuelMass", "fuel_mass_input")
    )
    dry_mass = streamlit.number_input(
        "launch angle (deg)",
        value=np.radians(2),
        format="%.5f",
        key="launch_angle_input",
        on_change=update_variable,
        args=("launchAngle", "launch_angle_input")
    )
    dry_mass = streamlit.number_input(
        "launch direction (deg)",
        value=np.radians(-135),
        format="%.5f",
        key="launch_direction_input",
        on_change=update_variable,
        args=("launchDirection", "launch_direction_input")
    )
    dry_mass = streamlit.number_input(
        "Center of gravity dry",
        value=1.0463,
        format="%.5f",
        key="cg_dry_input",
        on_change=update_variable,
        args=("CG_dry", "cg_dry_input")
    )
    dry_mass = streamlit.number_input(
        "Center of gravity wet",
        value=1.0774,
        format="%.5f",
        key="cg_wet_input",
        on_change=update_variable,
        args=("CG_wet", "cg_wet_input")
    )
    dry_mass = streamlit.number_input(
        "initial temperature (C)",
        value=20.53583333,
        format="%.8f",
        key="initial_temperature_input",
        on_change=update_variable,
        args=("InitialTemperature", "initial_temperature_input")
    )
    dry_mass = streamlit.number_input(
        "initial altitude (m)",
        value=492,
        format="%.5f",
        key="initial_altitude_input",
        on_change=update_variable,
        args=("initialAltitude", "initial_altitude_input")
    )
    dry_mass = streamlit.number_input(
        "initial pressure (Pa)",
        value=101550,
        format="%.5f",
        key="intial_pressure_input",
        on_change=update_variable,
        args=("initialPressure", "initial_pressure_input")
    )
    




if "results" in streamlit.session_state:
    df = pandas.DataFrame(streamlit.session_state.results)

    fig = plots.trajectory_3d_figure(df, color_by="time")
    streamlit.plotly_chart(fig, width="stretch") # 3D GRAPH
    # streamlit.dataframe(df)
    streamlit.text("Velocity VS Time")
    streamlit.line_chart(df, x="time", y="vel_up")
    streamlit.text("Altitude VS Time")
    streamlit.line_chart(df, x="time", y="pos_up")
    streamlit.text("Trajectory 2D VS Time")
    streamlit.scatter_chart(df, x="planar_distance", y="pos_up", size=5)
    streamlit.text("Vertical Force VS Time")
    streamlit.line_chart(df, x="time", y="force_up")

else:
    streamlit.info("Set your parameters and click **run simulation**.")
