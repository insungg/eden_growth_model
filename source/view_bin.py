import sys
from matplotlib import pyplot as plt
from decoder import read_bin

from draw_lattice import draw_lattice

def view_bin(file_path):
    params, landscapes, hotspots = read_bin(file_path)
    print("Parameters:")
    for key, value in params.items():
        print(f"{key}: {value}")

    while True:
        try:
            landscape_index = int(input(f"Enter the landscape index to plot (0 to {len(landscapes) - 1}, or -1 to exit): "))
            if landscape_index == -1:
                break
            if 0 <= landscape_index < len(landscapes):
                # draw_lattice_with_hotspots(landscapes[landscape_index], hotspots)
                draw_lattice(landscapes[landscape_index], hotspots)
                plt.show()  # Ensure plt.show() is called to display the plot
            else:
                print(f"Invalid index. Please enter a number between 0 and {len(landscapes) - 1}.")
        except ValueError:
            print("Invalid input. Please enter a valid number.")

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python view_bin.py <file_path>")
        sys.exit(1)

    file_path = sys.argv[1]
    view_bin(file_path)