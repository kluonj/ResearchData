import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

filename = "density_z.profile"

# Step 1: Parse the whole file into a list of (timestep, zlist, densitylist)
with open(filename) as f:
    lines = f.readlines()

frames = []
i = 0
while i < len(lines):
    line = lines[i].strip()
    if line.startswith('#') or not line:
        i += 1
        continue

    parts = line.split()
    if len(parts) == 3:
        try:
            timestep = int(parts[0])
            nchunks = int(parts[1])
            i += 1
            zlist = []
            densitylist = []
            for _ in range(nchunks):
                if i >= len(lines):
                    break
                dataline = lines[i].strip()
                datacols = dataline.split()
                if len(datacols) >= 4:
                    zlist.append(float(datacols[1]))
                    densitylist.append(float(datacols[3]))
                i += 1
            frames.append((timestep, zlist, densitylist))
        except ValueError:
            i += 1
    else:
        i += 1

# Step 2: Setup plot
fig, ax = plt.subplots(figsize=(6, 4))
line, = ax.plot([], [], 'o-', color='darkgreen')
ax.set_xlabel("z (Coord1)")
ax.set_ylabel("Density")
ax.set_title("Density Profile Over Time")
ax.grid(True)

# Create a text annotation to display timestep
timestep_text = ax.text(0.05, 0.95, '', transform=ax.transAxes, fontsize=12, verticalalignment='top')

def init():
    ax.set_xlim(min(frames[0][1]), max(frames[0][1]))
    ax.set_ylim(0, max(max(d) for _, _, d in frames))
    return line, timestep_text

def update(frame):
    timestep, zlist, densitylist = frame
    line.set_data(zlist, densitylist)
    timestep_text.set_text(f"Timestep: {timestep}")
    return line, timestep_text

# Step 3: Animate and save (optional)
ani = FuncAnimation(fig, update, frames=frames, init_func=init,
                    blit=True, interval=200)

# Optional: Save as MP4
# ani.save("density_profile.mp4", writer="ffmpeg", fps=5)
ani.save("density_profile.gif", writer="pillow", fps=2)

# Show the animation
#plt.show()

