import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.collections import PolyCollection
from matplotlib.colors import Normalize
from matplotlib import colormaps

def draw_lattice(lattice, hotspots, hex_radius=1, ax=None, fig_name=None):
    rows, cols = lattice.shape
    hex_width = np.sqrt(3) * hex_radius
    hex_height = 2 * hex_radius
    angles = np.linspace(np.pi / 6, 2 * np.pi + np.pi / 6, 7)
    hex_x = hex_radius * np.cos(angles)
    hex_y = hex_radius * np.sin(angles)
    
    # Generate hexagon coordinates
    x_indices, y_indices = np.meshgrid(np.arange(rows), np.arange(cols))
    center_x = x_indices * hex_width + (y_indices % 2) * (hex_width / 2)
    center_y = -(rows - 1 - y_indices) * (3 / 4 * hex_height) # grow from the bottom
    hex_centers = np.column_stack((center_x.flatten(), center_y.flatten()))
    hex_vertices = np.array([np.column_stack((hex_x + cx, hex_y + cy)) for cx, cy in hex_centers])
    # Colors for lattice
    color_map = {0: "white", 1: "red", 2: "yellow"}
    face_colors = np.array([color_map[lattice[i, j]] for i in range(rows) for j in range(cols)])
    hotspot_mask = hotspots.flatten() == 1
    edge_colors="none" if np.sqrt(rows*cols) >= 100 else "black"

    fig = None
    if ax is None:
        fig, ax = plt.subplots(figsize=(10, 10))
    # Regular hexagons
    ax.add_collection(PolyCollection(
        verts=hex_vertices,
        facecolors=face_colors,
        edgecolors=edge_colors
    ))
    # Hotspot hexagons with transparency
    ax.add_collection(PolyCollection(
        verts=hex_vertices[hotspot_mask],
        facecolors="black",
        edgecolors=edge_colors,
        alpha=0.4  # 50% transparency for hotspots
    ))
    ax.set_xlim(center_x.min() - hex_radius, center_x.max() + hex_radius)
    ax.set_ylim(center_y.min() - hex_radius, center_y.max() + hex_radius)
    ax.set_aspect('equal')
    ax.axis('off')
    if fig is not None: plt.show()
    if fig_name: fig.savefig(fig_name, dpi=300, bbox_inches="tight")

def draw_lattice_heatmap(density, hotspots=None, hex_radius=1, ax=None, color_map="turbo", vmin=0, vmax=0.5, fig_name=None):
    density = np.clip(density, vmin, vmax) # Probability > v_max is clipped to v_max 
    rows, cols = density.shape
    hex_width = np.sqrt(3) * hex_radius
    hex_height = 2 * hex_radius
    angles = np.linspace(np.pi / 6, 2 * np.pi + np.pi / 6, 7)
    hex_x = hex_radius * np.cos(angles)
    hex_y = hex_radius * np.sin(angles)
    
    # Generate hexagon coordinates
    x_indices, y_indices = np.meshgrid(np.arange(rows), np.arange(cols))
    center_x = x_indices * hex_width + (y_indices % 2) * (hex_width / 2)
    center_y = -(rows - 1 - y_indices) * (3 / 4 * hex_height) # grow from the bottom
    hex_centers = np.column_stack((center_x.flatten(), center_y.flatten()))
    hex_vertices = np.array([np.column_stack((hex_x + cx, hex_y + cy)) for cx, cy in hex_centers])
    
    # Normalize density values for color mapping
    norm = Normalize(vmin, vmax)
    cmap = colormaps[color_map]
    face_colors = cmap(norm(density.flatten()))
    
    fig = None
    if ax is None:
        fig, ax = plt.subplots(figsize=(10, 10))
    
    # Add hexagon patches
    ax.add_collection(PolyCollection(
        verts=hex_vertices,
        facecolors=face_colors,
        edgecolors="none"
    ))

    # Add hotspot hexagons with transparency
    if hotspots is not None:
        hotspot_mask = hotspots.flatten() == 1
        edge_colors="none" if np.sqrt(rows*cols) >= 100 else "black"
        ax.add_collection(PolyCollection(
        verts=hex_vertices[hotspot_mask],
        facecolors="black",
        edgecolors=edge_colors,
        alpha=0.4  # 50% transparency for hotspots
    ))
    
    ax.set_xlim(center_x.min() - hex_radius, center_x.max() + hex_radius)
    ax.set_ylim(center_y.min() - hex_radius, center_y.max() + hex_radius)
    ax.set_aspect('equal')
    ax.axis('off')
    
    if fig is not None: plt.show()
    if fig_name: fig.savefig(fig_name, dpi=300, bbox_inches="tight")