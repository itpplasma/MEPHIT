#!/usr/bin/env python3
"""
Plot large triangulation results from Fortran Delaunay triangulation implementation
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from pathlib import Path

def read_freefem_mesh(filename):
    """Read FreeFEM mesh format"""
    if not Path(filename).exists():
        print(f"Warning: {filename} not found, skipping")
        return None, None
    
    try:
        with open(filename, 'r') as f:
            lines = f.readlines()
        
        # Parse header
        header = lines[0].strip().split()
        npoints = int(header[0])
        ntriangles = int(header[1])
        
        if npoints == 0 or ntriangles == 0:
            print(f"Warning: {filename} has no valid data (points={npoints}, triangles={ntriangles})")
            return None, None
        
        # Parse points
        points = []
        for i in range(1, npoints + 1):
            parts = lines[i].strip().split()
            x, y = float(parts[0]), float(parts[1])
            points.append([x, y])
        
        # Parse triangles
        triangles = []
        for i in range(npoints + 1, npoints + 1 + ntriangles):
            parts = lines[i].strip().split()
            # Convert to 0-based indexing
            v1, v2, v3 = int(parts[0]) - 1, int(parts[1]) - 1, int(parts[2]) - 1
            triangles.append([v1, v2, v3])
        
        return np.array(points), np.array(triangles)
    
    except Exception as e:
        print(f"Error reading {filename}: {e}")
        return None, None

def plot_triangulation(points, triangles, title, ax, show_points=True, show_labels=False):
    """Plot a triangulation"""
    if points is None or triangles is None:
        ax.text(0.5, 0.5, 'No data', ha='center', va='center', transform=ax.transAxes)
        ax.set_title(title)
        return
    
    # Plot triangles
    for triangle in triangles:
        if len(triangle) == 3 and all(0 <= idx < len(points) for idx in triangle):
            tri_points = points[triangle]
            tri_patch = patches.Polygon(tri_points, fill=False, edgecolor='blue', linewidth=0.8, alpha=0.7)
            ax.add_patch(tri_patch)
    
    # Plot points
    if show_points:
        ax.scatter(points[:, 0], points[:, 1], color='red', s=30, zorder=5, alpha=0.8)
    
    # Add point labels for smaller cases
    if show_labels and len(points) <= 20:
        for i, (x, y) in enumerate(points):
            ax.annotate(f'{i+1}', (x, y), xytext=(5, 5), textcoords='offset points', 
                       fontsize=8, alpha=0.7)
    
    ax.set_xlim(points[:, 0].min() - 0.1, points[:, 0].max() + 0.1)
    ax.set_ylim(points[:, 1].min() - 0.1, points[:, 1].max() + 0.1)
    ax.set_aspect('equal')
    ax.set_title(title)
    ax.grid(True, alpha=0.3)

def main():
    """Generate plots for large triangulations"""
    
    # Test cases to plot
    test_cases = [
        ('large_random_10.msh', '10 Random Points'),
        ('large_random_25.msh', '25 Random Points'), 
        ('large_random_50.msh', '50 Random Points'),
        ('large_random_100.msh', '100 Random Points'),
        ('test_direct_fortran.msh', 'Simple Square'),
        ('test_complex_geometry.msh', 'Complex Hexagon')
    ]
    
    # Create figure with subplots
    fig, axes = plt.subplots(2, 3, figsize=(18, 12))
    axes = axes.flatten()
    
    for i, (filename, title) in enumerate(test_cases):
        points, triangles = read_freefem_mesh(filename)
        
        if points is not None and triangles is not None:
            npoints = len(points)
            ntriangles = len(triangles)
            title_with_stats = f"{title}\n{npoints} points, {ntriangles} triangles"
            
            # Show point labels only for smaller cases
            show_labels = npoints <= 10
            plot_triangulation(points, triangles, title_with_stats, axes[i], 
                             show_points=True, show_labels=show_labels)
        else:
            axes[i].text(0.5, 0.5, f'{title}\nNo data available', 
                        ha='center', va='center', transform=axes[i].transAxes)
            axes[i].set_title(title)
    
    plt.tight_layout()
    plt.savefig('plots/large_triangulations_overview.png', dpi=150, bbox_inches='tight')
    plt.savefig('plots/large_triangulations_overview.pdf', bbox_inches='tight')
    print("Generated plots/large_triangulations_overview.png")
    print("Generated plots/large_triangulations_overview.pdf")
    
    # Create individual detailed plots for the largest cases
    detailed_cases = [
        ('large_random_50.msh', '50 Random Points'),
        ('large_random_100.msh', '100 Random Points')
    ]
    
    for filename, title in detailed_cases:
        points, triangles = read_freefem_mesh(filename)
        
        if points is not None and triangles is not None:
            fig, ax = plt.subplots(1, 1, figsize=(12, 10))
            
            npoints = len(points)
            ntriangles = len(triangles)
            title_with_stats = f"{title}\n{npoints} points, {ntriangles} triangles"
            
            plot_triangulation(points, triangles, title_with_stats, ax, 
                             show_points=True, show_labels=False)
            
            # Add mesh quality info
            ax.text(0.02, 0.98, f'Triangles: {ntriangles}\nPoints: {npoints}\nRatio: {ntriangles/npoints:.1f}', 
                   transform=ax.transAxes, verticalalignment='top',
                   bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))
            
            output_name = filename.replace('.msh', '_detailed')
            plt.savefig(f'plots/{output_name}.png', dpi=150, bbox_inches='tight')
            plt.savefig(f'plots/{output_name}.pdf', bbox_inches='tight')
            print(f"Generated plots/{output_name}.png")
            plt.close()
    
    # Create a comparison plot showing triangle density
    fig, ax = plt.subplots(1, 1, figsize=(10, 6))
    
    sizes = []
    triangle_counts = []
    labels = []
    
    for filename, title in test_cases:
        points, triangles = read_freefem_mesh(filename)
        if points is not None and triangles is not None:
            sizes.append(len(points))
            triangle_counts.append(len(triangles))
            labels.append(title.split()[0])  # First word only
    
    if sizes:
        ax.scatter(sizes, triangle_counts, s=100, alpha=0.7, color='blue')
        for i, label in enumerate(labels):
            ax.annotate(label, (sizes[i], triangle_counts[i]), 
                       xytext=(5, 5), textcoords='offset points')
        
        # Add theoretical line (2n - 5 for planar graphs)
        x_theory = np.linspace(min(sizes), max(sizes), 100)
        y_theory = 2 * x_theory - 5
        ax.plot(x_theory, y_theory, 'r--', alpha=0.5, label='2n - 5 (theoretical)')
        
        ax.set_xlabel('Number of Points')
        ax.set_ylabel('Number of Triangles')
        ax.set_title('Triangulation Scaling: Points vs Triangles')
        ax.grid(True, alpha=0.3)
        ax.legend()
    
    plt.savefig('plots/triangulation_scaling.png', dpi=150, bbox_inches='tight')
    print("Generated plots/triangulation_scaling.png")
    plt.close()

if __name__ == '__main__':
    # Ensure plots directory exists
    Path('plots').mkdir(exist_ok=True)
    main()