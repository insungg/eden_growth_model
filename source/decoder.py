import re
import os
import numpy as np

def parse_filename(filename):
    pattern = re.compile(r".*/\w+_N_(?P<numSnapshots>\d+)_rows_(?P<rows>\d+)_cols_(?P<cols>\d+)_Ra_(?P<Ra>[\d.]+)_Rb_(?P<Rb>[\d.]+)_phi_(?P<phi>[\d.]+)_mu_(?P<mu>[-\d.]+)_nu_(?P<nu>[\d.]+)_s_(?P<s>[\d.]+)\.bin")
    match = pattern.match(filename)
    if not match:
        raise ValueError("Filename does not match expected pattern")
    return {key: float(value) if '.' in value else int(value) for key, value in match.groupdict().items()}

def read_bin(filename):
    with open(filename, 'rb') as file:
        data = file.read()
    
    params = parse_filename(filename)
    num_snapshots = params['numSnapshots']
    rows = params['rows']
    cols = params['cols']

    # Calculate the size of one snapshot
    snapshot_size = rows * cols  # 1 byte per uint8_t

    snapshots = []
    for n in range(num_snapshots):
        start = n * snapshot_size
        end = start + snapshot_size
        lattice = np.frombuffer(data[start:end], dtype=np.uint8).reshape((rows, cols)).astype(np.int32)
        snapshots.append(lattice)
        del lattice

    # Read the hotspots (assuming they are the same for all snapshots and stored at the end)
    hotspots = np.frombuffer(data[-snapshot_size:], dtype=np.uint8).reshape((rows, cols)).astype(np.int32)

    return params, snapshots, hotspots

def get_landscape(data_path, prefix, verbose=False):
    results = {}

    for filename in os.listdir(data_path):
        if filename.startswith(prefix) and filename.endswith(".bin"):
            file_path = os.path.join(data_path, filename)

            print(f"Processing file: {filename}")

            params, snapshots, hotspots = read_bin(file_path) 

            nu = params['nu']
            s = params['s']
            phase_params = (nu, s)

            # Initialize storage if this is the first file for this (nu, s)
            if phase_params not in results:
                results[phase_params] = {"snapshots": [], "hotspots": []}

            results[phase_params]["snapshots"].extend(snapshots)
            results[phase_params]["hotspots"].append(hotspots)  # Hotspots are the same for all snapshots in the file

            del snapshots, hotspots

    if verbose:
        # Display structure of results array
        print("Structure of results (sorted by phase_params):")
        sorted_phase_params = sorted(results.keys(), key=lambda x: (x[0], x[1]))
        for phase_params in sorted_phase_params:
            nu, s = phase_params
            data = results[phase_params]
            snapshots = data["snapshots"]
            hotspots = data["hotspots"]

            print(f"Phase Params: (nu={nu}, s={s})")
            print(f"  Number of Snapshots: {len(snapshots)}")
            if snapshots:
                # Print dimensions of the first snapshot to show structure
                print(f"  Snapshot Dimensions: {snapshots[0].shape}")
            else:
                print("  No snalshots available")

            print(f"  Number of Hotspots: {len(hotspots)}")
            if hotspots:
                # Print dimensions of the first hotspot to show structure
                print(f"  Hotspot Dimensions: {hotspots[0].shape}")
            else:
                print("  No hotspots available")

    return results