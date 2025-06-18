"""
@author: enmanuelduarte

Interactive GUI for discovering the behavior of a quantum particle in a two-dimensional infinite potential box. 
The box dimensions and quantum numbers can be changed for deeper insight. 
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import tkinter as tk
from tkinter import ttk

# Constants
hbar = 1.0545718e-34  # Reduced Planck's constant
m = 9.10938356e-31    # Mass of the electron

# Function to calculate the wavefunction
def wavefunction(nx, ny, Lx, Ly, x, y):
    return np.sin(nx * np.pi * x / Lx) * np.sin(ny * np.pi * y / Ly)

# Function to update the plot
def update_plot():
    nx = nx_slider.get()
    ny = ny_slider.get()
    Lx = float(Lx_entry.get())
    Ly = float(Ly_entry.get())
    plot_probability = prob_var.get()

    x = np.linspace(0, Lx, 100)
    y = np.linspace(0, Ly, 100)
    X, Y = np.meshgrid(x, y)

    if plot_probability:
        Z = (wavefunction(nx, ny, Lx, Ly, X, Y))**2
    else:
        Z = wavefunction(nx, ny, Lx, Ly, X, Y)

    ax.clear()
    ax.plot_surface(X, Y, Z, cmap='viridis')
    ax.set_title(f'Wavefunction for nx={nx}, ny={ny}')
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Wavefunction')

    canvas.draw()

# Create the main window
root = tk.Tk()
root.title("Quantum Particle in a 2D Box")

# Create the figure and axis for the plot
fig = plt.Figure(figsize=(6, 5), dpi=100)
ax = fig.add_subplot(111, projection='3d')

# Create the canvas for the plot
canvas = FigureCanvasTkAgg(fig, master=root)
canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=True)

# Create sliders for quantum numbers
nx_slider = tk.Scale(root, from_=1, to=10, orient=tk.HORIZONTAL, label='nx', command=lambda e: update_plot())
nx_slider.pack(fill=tk.X)
ny_slider = tk.Scale(root, from_=1, to=10, orient=tk.HORIZONTAL, label='ny', command=lambda e: update_plot())
ny_slider.pack(fill=tk.X)

# Create entry boxes for Lx and Ly
Lx_label = tk.Label(root, text="Lx:")
Lx_label.pack()
Lx_entry = tk.Entry(root)
Lx_entry.insert(0, "1.0")
Lx_entry.pack()

Ly_label = tk.Label(root, text="Ly:")
Ly_label.pack()
Ly_entry = tk.Entry(root)
Ly_entry.insert(0, "1.0")
Ly_entry.pack()

# Create a checkbox for plotting probability
prob_var = tk.BooleanVar()
prob_check = tk.Checkbutton(root, text="Plot Probability", variable=prob_var, command=update_plot)
prob_check.pack()

# Initialize the plot
update_plot()

# Start the main loop
root.mainloop()
