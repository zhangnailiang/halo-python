import matplotlib.pyplot as plt
import numpy as np
import spiceypy as sp
from process import rotational, Model
import json

fig = plt.figure(figsize=(10, 8))
ax = plt.axes(projection='3d')

### User's choice
MoonSphere = 0  # 1 if the Moon is drawn as a sphere, 0 for a point
RotationalF = 0  # Plot the sequential in the Earth-Moon rotational frame (0=inertial)
Converged = 0  # Plot the initial and converged trajectory (for optimization mode)
earth = 0  # Plot the earth (only in rotat. frame)
model = 0  # Plot the CR3BP trajectory (for optimization mode)
path = "input/statesNRHOCapstone.csv"

# Load the computed data
with open('output/ORBdata.json', 'r') as f:
    mat = json.load(f)

for e in list(mat["seq"].keys())[1:]:
    Ssat = np.array(mat["seq"][e]["XJ2000"])  # Convert to NumPy array
    T = np.array(mat["seq"][e]["t"])          # Convert to NumPy array
    if Converged:
        SConv = np.array(mat["seq"][e]["XC"])  # Convert to NumPy array if needed
    else:
        SConv = None

    # Process the change of frame and the model if needed
    if RotationalF:
        Ssat, SEarth, SConv = rotational(Ssat, T, Converged, SConv)
    if model:
        ModP = Model(path, RotationalF, T)
    Plot = Ssat

    # Plot the propagation
    ax.plot(Plot[:, 0], Plot[:, 1], Plot[:, 2], label="seq " + e, zorder=10)
    ax.scatter(Plot[0, 0], Plot[0, 1], Plot[0, 2])

    # Plot depending on user's choice
    if Converged:
        ax.plot(SConv[:, 0], SConv[:, 1], SConv[:, 2], "orange", zorder=10, label="Traj. converged")
        ax.scatter(SConv[0, 0], SConv[0, 1], SConv[0, 2], c="orange")
    if model:
        ax.plot(ModP[0], ModP[1], ModP[2], c="g", label="CR3BP")
        ax.scatter(ModP[0][0], ModP[1][0], ModP[2][0], c="g")
    if earth:
        if not RotationalF:
            sp.furnsh("input/de430.bsp")
            SEarth = np.zeros((len(T), 3))
            for i in range(len(T)):
                SEarth[i] = sp.spkezr("EARTH", T[i], "J2000", "NONE", "MOON")[0][:3]
        else:
            ax.scatter(-389703, 0, 0, c="b", label="Earth_CR3BP")
        ax.plot(SEarth[:, 0], SEarth[:, 1], SEarth[:, 2], c="c")
        ax.scatter(SEarth[0, 0], SEarth[0, 1], SEarth[0, 2], c="c")
if earth:
    ax.plot(SEarth[:, 0], SEarth[:, 1], SEarth[:, 2], c="c", label="Earth")

# Graph settings
ax.set_xlabel('X (km)')
ax.set_ylabel('Y (km)')
ax.set_zlabel('Z (km)')
if MoonSphere:  # Plot moon
    rM = 1738
    u = np.linspace(0, 2 * np.pi, 50)
    v = np.linspace(0, np.pi, 50)
    x = rM * np.outer(np.cos(u), np.sin(v))
    y = rM * np.outer(np.sin(u), np.sin(v))
    z = rM * np.outer(np.ones(np.size(u)), np.cos(v))
    ax.plot_surface(x, y, z, cmap="gray", zorder=0)
else:
    ax.scatter(0, 0, 0, c="gray", label="Moon")
plt.legend()

# Manually set equal axis limits
max_range = np.ptp([ax.get_xlim(), ax.get_ylim(), ax.get_zlim()]) / 2.0
mid_x = (ax.get_xlim()[0] + ax.get_xlim()[1]) / 2.0
mid_y = (ax.get_ylim()[0] + ax.get_ylim()[1]) / 2.0
mid_z = (ax.get_zlim()[0] + ax.get_zlim()[1]) / 2.0

ax.set_xlim(mid_x - max_range, mid_x + max_range)
ax.set_ylim(mid_y - max_range, mid_y + max_range)
ax.set_zlim(mid_z - max_range, mid_z + max_range)

sp.kclear()
plt.show()
