import io #agar

import streamlit
import pandas
import numpy as np

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
    # Keep results in memory (session_state), not a shared file on disk.
    streamlit.session_state.results = rocket.get_results()


def update_variable(variable_name, state_key):
    setattr(streamlit.session_state.sim, variable_name, streamlit.session_state[state_key])


with streamlit.sidebar:
    streamlit.button("run simulation", on_click=run_button_pressed)

    uploaded_file = streamlit.file_uploader("Choose a motor file", type="csv")
    if uploaded_file is not None:
        # Parse straight from the uploaded bytes - never touches disk, so
        # concurrent users can't clobber each other's motor file.
        text = uploaded_file.getvalue().decode("utf-8")
        streamlit.session_state.sim.set_motor_data(io.StringIO(text))

    dry_mass = streamlit.number_input(
        "dry mass (kg)",
        value=0.0605,
        format="%.5f",
        key="dry_mass_input",
        on_change=update_variable,
        args=("dryMass", "dry_mass_input")
    )
    fuel_mass = streamlit.number_input(
        "fuel mass (kg)",
        value=0.011,
        format="%.5f",
        key="fuel_mass_input",
        on_change=update_variable,
        args=("fuelMass", "fuel_mass_input")
    )
    radius = streamlit.number_input(
        "radius (m)",
        value=0.0127,
        format="%.5f",
        key="radius_input",
        on_change=update_variable,
        args=("radius", "radius_input")
    )
    drag_coefficient = streamlit.number_input(
        "drag coefficient",
        value=0.634,
        format="%.5f",
        key="drag_coefficient_input",
        on_change=update_variable,
        args=("dragCoefficient", "drag_coefficient_input")
    )
    parachute_area = streamlit.number_input(
        "parachute area (m^2)",
        value=0.0707,
        format="%.5f",
        key="parachute_area_input",
        on_change=update_variable,
        args=("parachuteArea", "parachute_area_input")
    )
    parachute_drag_coefficient = streamlit.number_input(
        "parachute drag coefficient",
        value=0.8,
        format="%.5f",
        key="parachute_drag_coefficient_input",
        on_change=update_variable,
        args=("parachuteDragCoefficient", "parachute_drag_coefficient_input")
    )
    rod_length = streamlit.number_input(
        "rod length (m)",
        value=1.0,
        format="%.5f",
        key="rod_length_input",
        on_change=update_variable,
        args=("rodLength", "rod_length_input")
    )
    launch_angle = streamlit.number_input(
        "launch angle (deg)",
        value=5.0,
        key="launch_angle_input",
        on_change=lambda: setattr(
            streamlit.session_state.sim, "launchAngle",
            np.radians(streamlit.session_state["launch_angle_input"])
        )
    )
    wind_speed = streamlit.number_input(
        "wind speed (m/s)",
        value=0.0,
        format="%.5f",
        key="wind_speed_input",
        on_change=update_variable,
        args=("windSpeed", "wind_speed_input")
    )
    rocket_length = streamlit.number_input(
        "rocket length (m)",
        value=0.425,
        format="%.5f",
        key="rocket_length_input",
        on_change=update_variable,
        args=("rocketLength", "rocket_length_input")
    )
    CG_dry = streamlit.number_input(
        "CG dry (m)",
        value=0.24,
        format="%.5f",
        key="cg_dry_input",
        on_change=update_variable,
        args=("CG_dry", "cg_dry_input")
    )
    CG_wet = streamlit.number_input(
        "CG wet (m)",
        value=0.26,
        format="%.5f",
        key="cg_wet_input",
        on_change=update_variable,
        args=("CG_wet", "cg_wet_input")
    )

if "results" in streamlit.session_state:
    df = pandas.DataFrame(streamlit.session_state.results)
    streamlit.line_chart(df, x="time", y="y_pos")
    streamlit.dataframe(df, use_container_width=True)
    streamlit.download_button(
        "Export as CSV",
        data=streamlit.session_state.sim.results_csv_bytes(),
        file_name="output.csv",
        mime="text/csv",
    )
else:
    streamlit.info("Set your parameters and click **run simulation**.")
