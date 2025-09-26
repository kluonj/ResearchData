import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.preprocessing import StandardScaler
from sklearn.mixture import GaussianMixture
from scipy.spatial import cKDTree


# ------------------- CONSTANTS & PATHS -------------------
dump_file     = "dump.desc"
output_dir    = "per_frame_results_filter"
os.makedirs(output_dir, exist_ok=True)

# Aluminium atomic mass ? grams per atom
mass_Al_g_per_atom = 26.9815 / 6.02214076e23  
# Å³ ? cm³
ANG3_TO_CM3 = 1e-24  

# ------------------- HELPERS -------------------
def parse_frame(lines):
    """
    Given the lines of one frame (from "ITEM: TIMESTEP" down to next),
    return ndarray of atom data and a header?column?index map.
    """
    header_map = {}
    records = []
    for i, line in enumerate(lines):
        if line.startswith("ITEM: ATOMS"):
            cols = line.split()[2:]
            header_map = {name: idx for idx, name in enumerate(cols)}
            for row in lines[i+1:]:
                if row.startswith("ITEM:"): break
                parts = row.split()
                rec = [int(parts[0]), int(parts[1])] + list(map(float, parts[2:]))
                records.append(rec)
            break
    return np.array(records), header_map

def write_dataframe_to_dump(df, reference_dumpfile, output_dumpfile):
    header_lines = []

    # Step 1: Read header from reference dump file (first snapshot only)
    with open(reference_dumpfile, 'r') as f:
        for line in f:
            header_lines.append(line)
            if line.strip().startswith("ITEM: ATOMS"):
                break  # Stop after reading the header

    # Step 2: Write header + new DataFrame data
    with open(output_dumpfile, 'w') as out:
        for line in header_lines:
            out.write(line)

        # Write atom lines using the same column order as in header
        col_names = header_lines[-1].strip().split()[2:]  # skip "ITEM: ATOMS"

        for _, row in df[col_names].iterrows():
            out.write(' '.join(str(v) for v in row.values) + '\n')

def write_dataframe_with_custom_columns(df, reference_dumpfile, output_dumpfile):
    header_lines = []

    # Step 1: Read header from reference dump file
    with open(reference_dumpfile, 'r') as f:
        for line in f:
            if line.strip().startswith("ITEM: ATOMS"):
                break  # Stop before writing ATOMS line
            header_lines.append(line)
    
    # Step 2: Start writing to new file
    with open(output_dumpfile, 'w') as out:
        # Write header lines before "ITEM: ATOMS"
        for line in header_lines:
            out.write(line)

        # Write new ITEM: ATOMS line based on DataFrame columns
        out.write("ITEM: ATOMS " + ' '.join(df.columns) + '\n')

        # Write each row of the DataFrame
        for _, row in df.iterrows():
            out.write(' '.join(str(v) for v in row.values) + '\n')

def extract_lammps_header(input_file, output_file):
    with open(input_file, 'r') as f_in, open(output_file, 'w') as f_out:
        for line in f_in:
            f_out.write(line)
            if line.strip().startswith("ITEM: ATOMS"):
                break

def majority_vote_refine(phases, positions, rcut=4.0, max_iter=3):
    """
    Iteratively refine labels using neighborhood majority voting.
    Parameters:
        phases: array of phase strings ['liquid', 'vapor', ...]
        positions: (N, 3) positions
        rcut: cutoff distance in Angstroms
        max_iter: number of iterations to smooth
    Returns:
        refined_phases: array of phase strings
    """
    phase_id_map = {'liquid': 0, 'interfacial': 1, 'vapor': 2}
    reverse_map = {v: k for k, v in phase_id_map.items()}
    label_ids = np.array([phase_id_map[p] for p in phases])

    tree = cKDTree(positions)
    refined = label_ids.copy()

    for _ in range(max_iter):
        new_labels = refined.copy()
        for i, pos in enumerate(positions):
            neighbors = tree.query_ball_point(pos, rcut)
            if len(neighbors) == 0:
                continue
            neighbor_labels = refined[neighbors]
            majority = np.bincount(neighbor_labels, minlength=3).argmax()
            new_labels[i] = majority
        refined = new_labels

    return np.array([reverse_map[l] for l in refined])



def refine_phase_labels(positions, initial_phases, rcut=4.0):
    """
    Refine phase labels based on neighbor majority voting.
    Arguments:
        positions: (N, 3) numpy array of atom positions
        initial_phases: array of strings ("liquid", "vapor", "interfacial")
        rcut: neighbor search cutoff in Å
    Returns:
        refined_phases: array of strings, same shape as input
    """
    tree = cKDTree(positions)
    refined_phases = initial_phases.copy()

    # Map for easy counting
    phase_id_map = {"liquid": 0, "interfacial": 1, "vapor": 2}
    #phase_id_map = {
    #        "vapor":0,
    #        "interfacial":1,
    #        "liquid":2
    #}

    reverse_map = {v: k for k, v in phase_id_map.items()}
    phase_ids = np.array([phase_id_map[p] for p in initial_phases])

    # First pass: liquid->interfacial if neighbors are mostly interfacial
    for i, pos in enumerate(positions):
        if phase_ids[i] == 0:   #liquid
            neighbors = tree.query_ball_point(pos, rcut)
            neighbor_phases = phase_ids[neighbors]
            counts = np.bincount(neighbor_phases, minlength=3)
            total = counts[0]  + counts[2]+ counts[1]
            if (counts[0] /(counts[0] + counts[2]+ counts[1]) < 0.2):  # mostly liquid
                refined_phases[i] = "vapor"


    phase_ids = np.array([phase_id_map[p] for p in refined_phases])
    # first pass: interfacial to vapor if neighbors are mostly  not interfacial
    for i, pos in enumerate(positions):
        if phase_ids[i] == 1:  # interfacial
            neighbors = tree.query_ball_point(pos, rcut)
            neighbor_phases = phase_ids[neighbors]
            counts = np.bincount(neighbor_phases, minlength=3)

            total = counts[0] + counts[2]+ counts[1]
            if (counts[0] / total  < 0.1):  #  if liquid atom around interfacial atom is really low, change it to vapor
                refined_phases[i] = "vapor"



    return np.array(refined_phases)


# ------------------- MAIN LOOP -------------------
with open(dump_file) as f:
    all_lines = f.readlines()

# find indices of each frame
frame_starts = [i for i, l in enumerate(all_lines) if l.strip()=="ITEM: TIMESTEP"]

timesteps   = []
liquid_rho  = []
vapor_rho   = []
liquid_pzz  = []
vapor_pzz   = []

for idx, start in enumerate(frame_starts):
    end = frame_starts[idx+1] if idx+1 < len(frame_starts) else len(all_lines)
    frame_lines = all_lines[start:end]

    # read timestep
    timestep = int(frame_lines[1].strip())
    if timestep < 40000:
        continue
    timesteps.append(timestep)

    # parse atom block
    data, hmap = parse_frame(frame_lines)
    ids  = data[:, hmap["id"]].astype(int)
    typeid = data[:, hmap["type"]].astype(int)
    x    = data[:, hmap["x"]].astype(float)
    y    = data[:, hmap["y"]].astype(float)
    z    = data[:, hmap["z"]].astype(float)
    CN   = data[:, hmap["c_cna"]].reshape(-1,1)
    Vor  = data[:, hmap["c_vor[1]"]].reshape(-1,1)
    pxx  = data[:, hmap["c_p[1]"]].astype(float)
    pyy  = data[:, hmap["c_p[2]"]].astype(float)
    pzz  = data[:, hmap["c_p[3]"]].astype(float)

    # stack & standardize
    X  = np.hstack([CN, Vor])
    Xs = StandardScaler().fit_transform(X)

    # 3-component GMM
    gmm    = GaussianMixture(n_components=3, tol=0.001, covariance_type="full", random_state=0)
    gmm.fit(Xs)
    labels = gmm.predict(Xs)
    means  = gmm.means_[:,0]        # mean CN for each component
    order  = np.argsort(means)     # [vapor_idx, interf_idx, liquid_idx]
    phase_map = {
        order[0]: "vapor",
        order[1]: "interfacial",
        order[2]: "liquid"
    }

    #mfilter=False
    positions = np.vstack([x, y, z]).T
    initial_phases = np.array([phase_map[l] for l in labels])
    #if mfilter==True:
    phases = refine_phase_labels(positions, initial_phases, rcut=4.0)
    #else:
    #    phases = initial_phases
    #phases = np.array([phase_map[l] for l in labels])
    #phases = majority_vote_refine(initial_phases, positions, rcut=4.0, max_iter=2)


    # ??? 3-Phase Classification Plot with LaTeX ???
    color_map = {"liquid":"blue", "vapor":"red", "interfacial":"gray"}
    colors    = [color_map[p] for p in phases]

    plt.figure(figsize=(7,6))
    plt.scatter(CN, Vor, c=colors, alpha=0.6, edgecolor="k", linewidth=0.2)
    plt.xlabel(r"Coordination Number $C_N$")
    plt.ylabel(r"Voronoi Volume $V_{\mathrm{vor}}\ (\mathrm{\AA}^3)$")
    plt.title(r"3-Phase Classification using GMM: $(C_N,\,V_{\mathrm{vor}})$")
    for phase in ["liquid","interfacial","vapor"]:
        plt.scatter([], [], c=color_map[phase], label=phase.capitalize(), alpha=0.6)
    plt.legend(title=r"Phase Identification")
    plt.tight_layout()
    plt.savefig(f"{output_dir}/frame_{timestep}_classification.png")
    plt.close()
    # ---------------------------------------------------------???

    # compute densities (g/cm³) for liquid & vapor
    rho = {}
    for p in ("liquid","vapor"):
        mask = (phases == p)
        if not mask.any():
            rho[p] = 0.0
        else:
            total_mass   = mask.sum() * mass_Al_g_per_atom
            total_volume = Vor[mask].sum() * ANG3_TO_CM3
            rho[p]       = total_mass / total_volume
    liquid_rho.append(rho["liquid"])
    vapor_rho.append(rho["vapor"])

    # compute pressure (kbar) for liquid & vapor
    pressure = {}
    for p in ("liquid","vapor"):
        mask = (phases == p)
        if not mask.any():
            pressure[p] = 0.0
        else:
            total_stress = pzz[mask].sum()
            total_volume = Vor[mask].sum() 
            pressure[p]       = total_stress / total_volume
    liquid_pzz.append(pressure["liquid"])
    vapor_pzz.append(pressure["vapor"])


    # save per-frame CSV
    df = pd.DataFrame({
        "id": ids,
        "type": typeid,
        "x": x,
        "y": y,
        "z": z,
        "CN": CN.flatten(),
        "VoronoiVolume": Vor.flatten(),
        "phase": phases
    })
    #phase_map = {'liquid': 0, 'vapor': 1, 'interfacial': 2}
    phase_map = {'liquid': 0,  'interfacial': 1, 'vapor': 2}
    df['phase_id'] = df['phase'].map(phase_map).astype(int)
    df = df.drop(columns='phase')
    

    #df.to_csv(f"{output_dir}/frame_{timestep}.dump", index=False)
    #write_dataframe_to_dump(df, dump_file,f"{output_dir}/frame_{timestep}.dump") 
    write_dataframe_with_custom_columns(df, dump_file,f"{output_dir}/frame_{timestep}.dump") 

    print(rf"Frame {timestep}: $\rho$_liquid={rho['liquid']:.4f}, $\rho$_vapor={rho['vapor']:.4f}")

# ------------------- SAVE DENSITY VS TIMESTEP -------------------
out = np.column_stack([timesteps, liquid_rho, vapor_rho])
np.savetxt(f"{output_dir}/density_vs_timestep.txt", out,
           header="timestep liquid_rho vapor_rho", fmt="%d %.6f %.6f")

# ------------------- SAVE pressure VS TIMESTEP -------------------
out = np.column_stack([timesteps, liquid_pzz, vapor_pzz])
np.savetxt(f"{output_dir}/pzz_vs_timestep.txt", out,
           header="timestep liquid_pzz vapor_pzz", fmt="%d %.6f %.6f")


# ------------------- DENSITY PLOT with LaTeX -------------------
timesteps  = np.array(timesteps)
liquid_rho = np.array(liquid_rho)
vapor_rho  = np.array(vapor_rho)

mean_liq, std_liq = liquid_rho.mean(), liquid_rho.std()
mean_vap, std_vap = vapor_rho.mean(),   vapor_rho.std()

plt.figure(figsize=(8,5))
plt.errorbar(timesteps, liquid_rho, yerr=std_liq, fmt="o", color="blue", label=r"$\rho_{\rm liq}$")
plt.errorbar(timesteps, vapor_rho,   yerr=std_vap, fmt="x", color="red",  label=r"$\rho_{\rm vap}$")

plt.axhline(mean_liq, color="blue", linestyle="--",
            label=rf"$\overline{{\rho}}_{{\rm liq}} = {mean_liq:.4f}\,\mathrm{{g/cm^3}}$")
plt.axhline(mean_vap, color="red",  linestyle="--",
            label=rf"$\overline{{\rho}}_{{\rm vap}} = {mean_vap:.4f}\,\mathrm{{g/cm^3}}$")

plt.xlabel(r"Timestep")
plt.ylabel(r"Density $\rho\ (\mathrm{g/cm^3})$")
plt.title(r"Liquid-Vapor Density vs. Timestep")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.savefig(f"{output_dir}/density_vs_timestep.png")
#plt.show()



# ------------------- pzz PLOT with LaTeX -------------------
timesteps  = np.array(timesteps)
liquid_pzz = np.array(liquid_pzz)
vapor_pzz  = np.array(vapor_pzz)

mean_liq, std_liq = liquid_pzz.mean(), liquid_pzz.std()
mean_vap, std_vap = vapor_pzz.mean(),   vapor_pzz.std()

plt.figure(figsize=(8,5))
plt.errorbar(timesteps, liquid_pzz, yerr=std_liq, fmt="o", color="blue", label=r"$P_{zz}^{\rm liq}$")
plt.errorbar(timesteps, vapor_pzz,   yerr=std_vap, fmt="x", color="red",  label=r"$P_{zz}^{\rm vap}$")

#plt.axhline(mean_liq, color="blue", linestyle="--",
#            label=rf"$\overline{p_{zz}^{liq}} = {mean_liq:.4f}\,\mathrm{kbar}$")
#plt.axhline(mean_vap, color="red",  linestyle="--",
#            label=rf"$\overline{p_{zz}^{vap}} = {mean_vap:.4f}\,\mathrm{kbar}$")
plt.axhline(mean_liq, color="blue", linestyle="--",
            label=rf"$\overline{{p_{zz}}}_{{\rm liq}} = {mean_liq:.4f}\,\mathrm{{kbar}}$")
plt.axhline(mean_vap, color="red",  linestyle="--",
            label=rf"$\overline{{p_{zz}}}_{{\rm vap}} = {mean_vap:.4f}\,\mathrm{{kbar}}$")


plt.xlabel(r"Timestep")
plt.ylabel(r"$P_{zz}\ (\mathrm{g/cm^3})$")
plt.title(r"Liquid-Vapor pzz vs. Timestep")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.savefig(f"{output_dir}/pzz_vs_timestep.png")
#plt.show()


ex=np.array([mean_liq, std_liq, mean_vap, std_vap])
np.savetxt("pzz_filter.dat", np.transpose(ex))
