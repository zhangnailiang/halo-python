import json
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import spiceypy as sp
from tqdm import tqdm

### Load the computed data of propagation from JSON file
with open('output/ORBdata1.json', 'r') as file:
    data = json.load(file)

# Assuming the JSON structure contains keys similar to the .mat file
XJ2000 = np.array(data["seq"]["a"]["XJ2000"])
T = np.array(data["seq"]["a"]["t"])
XDF = pd.DataFrame(XJ2000, columns=["X", "Y", "Z", "VX", "VY", "VZ"])
start, stop = T[0], T[-1]

### Load Spice kernels
sp.furnsh("input/LRO_ES_90_202003_GRGM900C_L600.bsp")
sp.furnsh("input/naif0012.tls")  # Time management: leap seconds kernel

### User's choice
g = 0  # 1 to plot the trajectory
h = 0  # 1 if the data from the spacecraft comes from a Horizon file (JSON in this case)
SC = "LRO"  # Spacecraft (SC) if the Spice database is used

### Load data from a Horizon file if enabled
if h:
    with open("input/CAPSTONE_25_00__01_12.json", "r") as file:
        Hrz = pd.DataFrame(json.load(file))

### Plot trajectory with Moon
if g:
    fig = plt.figure()
    ax = plt.axes(projection='3d')

    # Plot Moon
    rM = 1738
    u = np.linspace(0, 2 * np.pi, 50)
    v = np.linspace(0, np.pi, 50)
    x = rM * np.outer(np.cos(u), np.sin(v))
    y = rM * np.outer(np.sin(u), np.sin(v))
    z = rM * np.outer(np.ones(np.size(u)), np.cos(v))
    ax.plot_surface(x, y, z, cmap="gray", zorder=1)

    # Plot trajectory
    ax.plot(XJ2000[:, 0], XJ2000[:, 1], XJ2000[:, 2], "r", zorder=10, label="Trajectory")
    ax.set_xlabel('X (km)')
    ax.set_ylabel('Y (km)')
    ax.set_zlabel('Z (km)')

    # Plot SC initial position
    state = XJ2000[0]
    ax.plot(state[0], state[1], state[2], "b.", zorder=10, label="Initial Pos")  # Starting point

    # Plot SC final position
    if h:
        state = Hrz.iloc[-1]
        state = [state["X"], state["Y"], state["Z"]]
    else:
        state = sp.spkezr(SC, stop, "J2000", "NONE", "MOON")[0]
    ax.plot(state[0], state[1], state[2], "c.", zorder=10, label="Final Pos")  # Ending point True

    # Plot final position of propagation
    state = XJ2000[-1]
    ax.plot(state[0], state[1], state[2], "y.", zorder=10, label="Final Pos Propagation")  # Ending point All
    ax.legend()
    plt.axis("equal")

### Acquire data for the errors plot
Tstep = T[1] - T[0]
Tspan = stop - start
Err = []
for _ in range(6):
    Err.append([])

for i in tqdm(range(len(XJ2000))):
    Time = T[i]
    # Plot SC final position
    if h:
        state = Hrz.iloc[i]
        state = [state["X"], state["Y"], state["Z"], state["VX"], state["VY"], state["VZ"]]
    else:
        state = sp.spkezr(SC, Time, "J2000", "NONE", "MOON")[0]
    for j in range(6):
        Err[j].append(state[j] - XJ2000[i, j])

### Plot position and velocity errors graph
fig2 = plt.figure()
axes = fig2.subplots(2, 1)
Error3D = [np.sqrt((np.array(Err[:3])**2).sum(0)) * 1e3, np.sqrt((np.array(Err[3:])**2).sum(0)) * 1e5]
for i in range(2):
    axes[i].plot((T - T[0]) / 3600, Error3D[i],
                 label="RMSE = " + '{:.2f}'.format(np.sqrt((Error3D[i]**2).sum() / len(Error3D[i]))) + [" m", " cm/s"][i])
    axes[i].set_ylabel(f'error in {["position (m)", "velocity (cm/s)"][i]}')
    axes[i].legend()
axes[1].set_xlabel('hours')

### Close the kernels used and show graphs
sp.kclear()
plt.show()
