#plots.py - 3D trajectory figure for RocketSim
#
# Builds a Plotly figure from RocketSim.get_results(). Works in Streamlit
# (st.plotly_chart) or on its own (python plots.py opens it in a browser).
#
# Axes: x = East, y = North, z = Up. Plotting North on x would mirror the
# view, so this keeps it looking like a map seen from the pad.

import numpy as np
import plotly.graph_objects as go

# Columns that make sense as a line colour, with their colour-bar titles
COLOR_OPTIONS = {
    "speed": "Speed (m/s)",
    "mach": "Mach",
    "pos_up": "Altitude (m)",
    "time": "Time (s)",
    "aoa": "Angle of attack (deg)",
}


def trajectory_3d_figure(results, color_by="speed", height=600):
    """3D flight path coloured by `color_by`, with a ground shadow and
    markers at the pad, apogee, and landing point.

    results  : dict of lists (RocketSim.get_results()) or a pandas DataFrame
    color_by : any column name; COLOR_OPTIONS lists the useful ones
    """
    east = np.asarray(results["pos_east"])
    north = np.asarray(results["pos_north"])
    up = np.asarray(results["pos_up"])
    t = np.asarray(results["time"])
    c = np.asarray(results[color_by])
    speed = np.asarray(results["speed"])

    apo = int(np.argmax(up))
    color_title = COLOR_OPTIONS.get(color_by, color_by)

    fig = go.Figure()

    # Shadow of the path on the ground - makes depth much easier to judge
    fig.add_trace(go.Scatter3d(
        x=east, y=north, z=np.zeros_like(up),
        mode="lines",
        line=dict(color="rgba(150,150,150,0.6)", width=2, dash="dash"),
        hoverinfo="skip",
        showlegend=False,
    ))

    # Flight path, coloured by the chosen column
    fig.add_trace(go.Scatter3d(
        x=east, y=north, z=up,
        mode="lines",
        line=dict(
            width=6,
            color=c,
            colorscale="Viridis",
            colorbar=dict(title=dict(text=color_title), thickness=12,
                          len=0.6),
        ),
        customdata=np.column_stack((t, speed, c)),
        hovertemplate=(
            "t = %{customdata[0]:.2f} s<br>"
            "East %{x:.1f} m<br>North %{y:.1f} m<br>Alt %{z:.1f} m<br>"
            "Speed %{customdata[1]:.1f} m/s<br>"
            + color_title + ": %{customdata[2]:.3g}<extra></extra>"
        ),
        showlegend=False,
    ))

    # Pad, apogee, landing
    fig.add_trace(go.Scatter3d(
        x=[0, east[apo], east[-1]],
        y=[0, north[apo], north[-1]],
        z=[0, up[apo], 0],
        mode="markers+text",
        marker=dict(size=5, color=["#888780", "#eb6834", "#1baf7a"]),
        text=["Pad", f"Apogee {up[apo]:.0f} m", "Landing"],
        textposition="top center",
        hovertemplate="%{text}<br>East %{x:.1f} m<br>North %{y:.1f} m"
                      "<br>Alt %{z:.1f} m<extra></extra>",
        showlegend=False,
    ))

    # Give East and North the same span (the larger of the two, padded), so a
    # flight with no East motion still gets a readable East axis instead of
    # one squashed to zero width. Metres stay equal on all three axes.
    x_range, y_range, z_range = _equal_scale_ranges(east, north, up)
    span = x_range[1] - x_range[0]
    fig.update_layout(
        height=height,
        margin=dict(l=0, r=0, t=0, b=0),
        scene=dict(
            xaxis=dict(title="East (m)", range=x_range),
            yaxis=dict(title="North (m)", range=y_range),
            zaxis=dict(title="Altitude (m)", range=z_range),
            aspectmode="manual",
            aspectratio=dict(x=1, y=1, z=(z_range[1] - z_range[0]) / span),
            camera=dict(eye=dict(x=-1.4, y=-1.5, z=0.8)),
        ),
    )
    return fig


def _equal_scale_ranges(east, north, up, pad=0.08, min_ground=0.5):
    """Axis ranges where East and North share one span (always including
    the pad at 0,0), and altitude runs from the ground to apogee. The ground
    is at least `min_ground` x apogee wide, so a near-vertical flight shows
    some empty ground instead of becoming a tall, thin sliver."""
    lo_e, hi_e = min(east.min(), 0.0), max(east.max(), 0.0)
    lo_n, hi_n = min(north.min(), 0.0), max(north.max(), 0.0)
    span = max(hi_e - lo_e, hi_n - lo_n, min_ground * up.max(), 10.0) * (
        1 + 2 * pad)
    mid_e, mid_n = (lo_e + hi_e) / 2, (lo_n + hi_n) / 2
    z_top = max(up.max(), 1.0) * (1 + pad)
    return ([mid_e - span / 2, mid_e + span / 2],
            [mid_n - span / 2, mid_n + span / 2],
            [0.0, z_top])


if __name__ == "__main__":
    # Standalone demo: run the default rocket with some wind and open the
    # figure in a browser. Run from the rocketsimweb folder.
    import math
    import sim

    s = sim.RocketSim(launchDirection=math.radians(60), windSpeed=3.0,
                      windSeed=1)
    s.run_simulation()
    trajectory_3d_figure(s.get_results()).show()
