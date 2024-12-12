import random
import os
import sys
import numpy as np

def generate_hotspots(rows, cols, Ra, Rb, phi, centers=None):
    lattice = np.zeros((rows, cols), dtype=int)
    if centers is None:
        num_hotspots = int(phi * rows * cols / (1 + 3 * Ra * (Rb + 1)))
        cx = np.random.randint(0, rows, size=num_hotspots)
        cy = np.random.randint(0, cols, size=num_hotspots)
        centers = list(zip(cx, cy))

    def to_axial(row, col):
        q = col - (row // 2)
        r = row
        return q, r

    def to_real(q, r):
        col = q + (r // 2)
        row = r
        return row, col

    for cx, cy in centers:
        center_q, center_r = to_axial(cx, cy)
        for dq in range(-Ra, Ra + 1):
            for dr in range(-Rb, Rb + 1):
                if (dq / Ra) ** 2 + (dr / Rb) ** 2 + ((dq + dr) / Ra) ** 2 <= 1:
                    row, col = to_real(center_q + dq, center_r + dr)
                    if 0 <= row < rows and 0 <= col < cols:
                        lattice[row, col] = 1
    return lattice

def make_parameters_file(file_name, mutant_fraction, output_dir):
    """Generate parameter files with specified mutant fraction and save to the output directory."""
    os.makedirs(output_dir, exist_ok=True)  # Ensure the output directory exists

    # Define parameters
    num_landscapes = 50
    num_snapshots = 500
    parallel_on = True
    num_cores = 10
    rows = 500
    cols = 500
    Ra = 5
    Rb = 5
    phi = 0.25
    mu = -1
    nu_values = [0, 2]
    s_values = [0]


    for i in range(1, num_landscapes + 1):
        # Generate initial condition for the first row with the specified mutant fraction
        num_mutants = int(cols * mutant_fraction)
        num_others = cols - num_mutants

        initial_condition = [2] * num_mutants + [1] * num_others
        random.shuffle(initial_condition)  # Randomize the placement of 1s and 2s
        # initial_condition = [1] * (cols // 2) + [2] + [1] * (cols - cols // 2 - 1) # Mutant in the middle

        # Format the initial condition as a comma-separated string
        formatted_condition = ','.join(map(str, initial_condition))

        # Format nu and s values as comma-separated strings
        formatted_nu_values = ','.join(map(str, nu_values))
        formatted_s_values = ','.join(map(str, s_values))

        # Generate initial hotspots
        # center = [[Rb, cols // 2]]
        # z = 180
        # centers = [[Rb + i * z - 1, cols // 2] for i in range(rows // z)]
        # print(centers)
        # initial_hotspots = generate_hotspots(rows, cols, Ra, Rb, phi)
        # formatted_hotspots = ';'.join(','.join(map(str, row)) for row in initial_hotspots)

        # Define the file content
        content = f"""numSnapshots={num_snapshots}
parallelOn={str(parallel_on).lower()}
numCores={num_cores}
fileName={file_name}hs{i}
rows={rows}
cols={cols}
Ra={Ra}
Rb={Rb}
phi={phi}
mu={mu}
nu={formatted_nu_values}
s={formatted_s_values}
initialCondition={formatted_condition}
"""

        # Write the content to a file
        output_file_path = os.path.join(output_dir, f"{file_name}hs{i}.txt")
        with open(output_file_path, "w") as file:
            file.write(content)

    print(f"Parameter files generated successfully in {output_dir}.")

if __name__ == "__main__":
    # Parse command-line arguments
    if len(sys.argv) != 4:
        print("Usage: python script.py <file_name> <mutant_fraction> <output_directory>")
        sys.exit(1)

    file_name = sys.argv[1]
    try:
        mutant_fraction = float(sys.argv[2])
    except ValueError:
        print("Error: mutant_fraction must be a floating-point number.")
        sys.exit(1)

    output_directory = sys.argv[3]
    make_parameters_file(file_name, mutant_fraction, output_directory)