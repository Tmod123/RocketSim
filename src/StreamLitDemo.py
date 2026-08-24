import streamlit as st
import plotly.express as px

st.title("Rocket Flight Model")

mass = st.slider(
    "Rocket Mass (g)",
    300,
    1500,
    650
)

cd = st.slider(
    "Drag Coefficient",
    0.2,
    1.0,
    0.5
)

# Run rocket model here
time, altitude, velocity = run_simulation(mass, cd)

fig = px.line(
    x=time,
    y=altitude,
    labels={
        "x": "Time (s)",
        "y": "Altitude (m)"
    }
)

st.plotly_chart(fig)

