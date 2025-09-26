import numpy as np
import matplotlib.pyplot as plt

filename = "density_z.profile"

with open(filename) as f:
    lines = f.readlines()

# Step 1: Find the last data block header
header_index = None
nchunks = None
timestep = None

for i in reversed(range(len(lines))):
    line = lines[i].strip()
    if not line or line.startswith('#'):
        continue
    parts = line.split()
    if len(parts) == 3:
        try:
            timestep = int(parts[0])
            nchunks = int(parts[1])
            header_index = i
            break
        except ValueError:
            continue

if header_index is None:
    raise RuntimeError("No valid data block header found.")

# Step 2: Read following nchunks lines
block = lines[header_index + 1 : header_index + 1 + nchunks]

zlist = []
densitylist = []

for line in block:
    parts = line.strip().split()
    if len(parts) >= 4:
        z = float(parts[1])  # Coord1
        density = float(parts[3])  # density/mass
        zlist.append(z)
        densitylist.append(density)

# Step 3: Plot
plt.figure(figsize=(6, 4))
plt.plot(zlist, densitylist, marker='o', linestyle='-', color='steelblue')
plt.xlabel("z (Coord1) [distance]")
plt.ylabel("Mass Density [density/mass]")
plt.title(f"Density Profile (Final Timestep: {timestep})")
plt.grid(True)
plt.tight_layout()
plt.savefig("density_z_profile.png", dpi=300)

#plt.show()

