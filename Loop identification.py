import os
from Bio.PDB import MMCIFParser, DSSP

# ==========================================================
# CONFIGURATION
# ==========================================================

cif_folder = "/home/fernando/Descargas/loop-protein-54"
output_file = "loop_results.txt"

# DSSP secondary structure codes considered as regular structure
# H = Alpha helix
# G = 3₁₀ helix
# I = Pi helix
# E = Beta strand
# B = Isolated beta bridge
#
# Any other DSSP state (blank, T, S, etc.) is considered part of a loop.

STRUCTURED_STATES = {"H", "G", "I", "E", "B"}

# ==========================================================
# FUNCTIONS
# ==========================================================

def get_longest_loop(structure, cif_file, model_id=0):
    """
    Identify the longest continuous loop region in a protein structure.

    Loop residues are defined as any DSSP assignment that is not
    H, G, I, E, or B.
    """

    model = structure[model_id]
    dssp = DSSP(model, cif_file, dssp="mkdssp")

    # DSSP returns tuples containing residue information and
    # secondary structure assignments.
    ss_data = [(res[0], res[1]) for res in dssp.property_list]

    loops = []
    current_loop = []

    for residue, ss in ss_data:

        if ss not in STRUCTURED_STATES:
            current_loop.append(residue)

        else:
            if current_loop:
                loops.append(current_loop)
                current_loop = []

    if current_loop:
        loops.append(current_loop)

    if not loops:
        return None

    longest_loop = max(loops, key=len)

    chain_id = longest_loop[0][0][0]
    start = longest_loop[0][1][1]
    end = longest_loop[-1][1][1]

    sequence = "".join(
        residue[0].get_resname()
        for residue in longest_loop
    )

    return chain_id, start, end, len(longest_loop), sequence


# ==========================================================
# MAIN ANALYSIS
# ==========================================================

parser = MMCIFParser(QUIET=True)
results = []

for filename in os.listdir(cif_folder):

    if filename.endswith(".cif"):

        cif_path = os.path.join(cif_folder, filename)

        try:

            structure = parser.get_structure(filename, cif_path)

            loop_info = get_longest_loop(structure, cif_path)

            if loop_info:

                chain, start, end, length, sequence = loop_info

                results.append(
                    f"{filename}\t{chain}\t{start}-{end}\t{length}\t{sequence}"
                )

            else:

                results.append(
                    f"{filename}\tNo loop region detected"
                )

        except Exception as e:

            results.append(
                f"{filename}\tError: {e}"
            )


# ==========================================================
# WRITE OUTPUT
# ==========================================================

with open(output_file, "w") as out:

    out.write("Structure\tChain\tRange\tLength\tSequence\n")

    for line in results:
        out.write(line + "\n")


print(
    f"Analysis completed successfully.\n"
    f"Results written to: {output_file}"
)
