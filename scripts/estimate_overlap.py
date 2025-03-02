import os
import glob

# Define input files
input_files = glob.glob("../repliseq_data/*lagging") + glob.glob("../repliseq_data/*leading") + glob.glob("../repliseq_data/*deciles.bed")

# Ensure output directories exist
os.makedirs("outputs", exist_ok=True)
os.makedirs("outputs_comb", exist_ok=True)

# First round of intersections (with -wb flag)
for file in input_files:
    output_file = f"outputs/{os.path.basename(file)}.i_wb.G4"
    command = f"bedtools intersect -a chm13v2.0_quadron_matches_corr_cordinates.bed.hg19.bed -b {file} -wb > {output_file}"
    os.system(command)
    print(command)

# Second round of intersections (with -u flag)
input_files = glob.glob("../repliseq_data/*lagging") + glob.glob("../repliseq_data/*leading")

for file in input_files:
    output_file = f"outputs_comb/{os.path.basename(file)}.i_u.G4"
    command = f"bedtools intersect -a chm13v2.0_quadron_matches_corr_cordinates.bed.hg19.bed -b {file} -u > {output_file}"
    os.system(command)

# Final intersection with deciles
for file in glob.glob("outputs_comb/*"):
    cell = os.path.basename(file).split(".i_u")[0].split("_")[0]
    decile_file = f"../{cell}_replitime.deciles.bed"

    if os.path.exists(decile_file):
        output_file = f"{file}.i_wb.{os.path.basename(decile_file)}"
        command = f"bedtools intersect -a {file} -b {decile_file} -wb > {output_file}"
        os.system(command)
    else:
        print(f"Warning: Decile file not found for cell type {cell}: {decile_file}")
