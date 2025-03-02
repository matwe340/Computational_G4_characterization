import os
import glob
import numpy as np
import matplotlib.pyplot as plt

# Define constants
OUTPUT_DIR = "out2"
ASSOCIATIONS_DIR = "associations"

os.makedirs(OUTPUT_DIR, exist_ok=True)
os.makedirs(ASSOCIATIONS_DIR, exist_ok=True)

def reader(path):
    """Reads a tab-separated file into a list of lists."""
    with open(path, "r") as file:
        return [line.strip().split('\t') for line in file.readlines()]


def process_leading_lagging():
    """Processes leading and lagging strand counts."""
    file_types = glob.glob("../*_RepliStrand.lagging")
    cell_lines = [os.path.basename(k).split("_")[0] for k in file_types]

    with open("Leading_Lagging_Asymmetry.txt", "w") as datafile:
        for cell in cell_lines:
            data_l1 = reader(os.path.join(OUTPUT_DIR, f"{cell}_RepliStrand.leading.i_wb.G4"))
            data_l2 = reader(os.path.join(OUTPUT_DIR, f"{cell}_RepliStrand.lagging.i_wb.G4"))

            leading, lagging = 0, 0
            for entry in data_l1 + data_l2:
                if entry[3] == "+":
                    leading += 1
                else:
                    lagging += 1

            print(cell, leading, lagging)
            datafile.write(f"{cell}\t{leading}\t{lagging}\n")


def process_deciles():
    """Processes replication timing deciles."""
    file_types = glob.glob("../*_replitime.deciles.bed")
    cell_lines = [os.path.basename(k).split("_")[0] for k in file_types]

    decile_dicts = {
        "All": {},
        "GGGA": {},
        "GGGT": {},
        "C9orf72": {},
        "SVA": {},
        "L1": {}
    }

    with open("Deciles.txt", "w") as datafile:
        for cell in cell_lines:
            data_l1 = reader(os.path.join(OUTPUT_DIR, f"{cell}_replitime.deciles.bed.i_wb.G4"))

            for dec in range(1, 11):
                found = sum(1 for k in data_l1 if int(k[-1]) == dec)
                datafile.write(f"{cell}\t{dec}\t{found}\n")

                # Update general decile dictionary
                decile_dicts["All"].setdefault(dec, []).append(found)

                # Update motif-specific counts
                for motif in ["GGGA", "GGGT", "C9orf72", "SVA", "L1"]:
                    found_motif = sum(1 for k in data_l1 if int(k[-1]) == dec and (k[4] == motif or k[4] == f"{motif}_rev_comp"))
                    decile_dicts[motif].setdefault(dec, []).append(found_motif)

    return decile_dicts


def plot_deciles(decile_dicts):
    """Generates and saves bar plots for deciles."""
    def plot_data(label, data_dict, filename):
        means = [np.mean(data_dict[v]) for v in range(1, 11)][::-1]
        stds = [np.std(data_dict[v]) for v in range(1, 11)][::-1]

        plt.bar(range(1, 11), means, yerr=stds, color="grey", ecolor="black")
        plt.grid()
        plt.ylabel("Occurrences", fontsize=16)
        plt.xlabel("Deciles", fontsize=16)
        plt.tight_layout()
        plt.savefig(filename)
        plt.close()

    # Plot for all motifs
    plot_data("All Motifs", decile_dicts["All"], "RepliSeq_All_Motifs_not_exact.png")

    # Density plot for all motifs
    total_np = sum([np.mean(decile_dicts["All"][v]) for v in range(1, 11)][::-1])
    density_means = [np.mean(decile_dicts["All"][v]) / total_np for v in range(1, 11)][::-1]
    density_stds = [np.std([mm / total_np for mm in decile_dicts["All"][v]]) for v in range(1, 11)][::-1]

    plt.bar(range(1, 11), density_means, yerr=density_stds, color="grey", ecolor="black")
    plt.grid()
    plt.ylabel("Density", fontsize=16)
    plt.xlabel("Deciles", fontsize=16)
    plt.tight_layout()
    plt.savefig("RepliSeq_All_Motifs_density_not_exact.png")
    plt.close()

    # Individual motif plots
    for motif in ["GGGA", "GGGT", "C9orf72", "SVA", "L1"]:
        plot_data(motif, decile_dicts[motif], f"RepliSeq_All_Motifs_{motif}_not_exact.png")

        total_motif = sum([np.mean(decile_dicts[motif][v]) for v in range(1, 11)][::-1])
        density_means = [np.mean(decile_dicts[motif][v]) / total_motif for v in range(1, 11)][::-1]
        density_stds = [np.std([mm / total_motif for mm in decile_dicts[motif][v]]) for v in range(1, 11)][::-1]

        plt.bar(range(1, 11), density_means, yerr=density_stds, color="grey", ecolor="black")
        plt.grid()
        plt.ylabel("Density", fontsize=16)
        plt.xlabel("Deciles", fontsize=16)
        plt.tight_layout()
        plt.savefig(f"RepliSeq_All_Motifs_{motif}_density_not_exact.png")
        plt.close()


def main():
    """Main function to execute all steps."""
    process_leading_lagging()
    decile_dicts = process_deciles()
    plot_deciles(decile_dicts)

if __name__ == "__main__":
    main()

