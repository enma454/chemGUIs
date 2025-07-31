#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 31 00:10:45 2025

@author: enmanuelduarte
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.animation import FuncAnimation
import tkinter as tk
from tkinter import ttk

class CoolingApp:
    def __init__(self, master):
        self.master = master
        self.master.title("Coffee Cooling Simulation (Newton's Law)")

        # Time array
        self.t = np.linspace(0, 60, 300)  # simulate 60 minutes

        # Initial parameters
        self.k = tk.DoubleVar(value=0.03)
        self.T_initial = tk.DoubleVar(value=90.0)
        self.T_ambient = tk.DoubleVar(value=20.0)

        # Set up plot
        self.fig, self.ax = plt.subplots(figsize=(6, 4))
        self.line, = self.ax.plot([], [], lw=2, color='brown')
        self.ax.set_title("Cooling of Coffee Over Time")
        self.ax.set_xlabel("Time (minutes)")
        self.ax.set_ylabel("Temperature (°C)")
        self.ax.set_ylim(0, 100)
        self.ax.grid()

        self.canvas = FigureCanvasTkAgg(self.fig, master=self.master)
        self.canvas.get_tk_widget().grid(row=0, column=0, columnspan=3)

        # Sliders
        self.create_slider("Cooling Rate (k)", self.k, 0.005, 0.1, 0)
        self.create_slider("Initial Temp (°C)", self.T_initial, 30, 100, 1)
        self.create_slider("Ambient Temp (°C)", self.T_ambient, 0, 40, 2)

        # Animation setup
        self.anim = FuncAnimation(self.fig, self.update_plot, frames=len(self.t),
                                  interval=30, blit=False, repeat=True)

    def create_slider(self, label, var, min_val, max_val, row):
        ttk.Label(self.master, text=label).grid(row=row+1, column=0, sticky='w')
        slider = ttk.Scale(self.master, variable=var, from_=min_val, to=max_val,
                           orient=tk.HORIZONTAL, command=lambda val: self.reset_plot())
        slider.grid(row=row+1, column=1, columnspan=2, sticky='we')
        self.master.grid_rowconfigure(row+1, pad=5)

    def cooling_function(self, t):
        T0 = self.T_initial.get()
        Ta = self.T_ambient.get()
        k = self.k.get()
        return Ta + (T0 - Ta) * np.exp(-k * t)

    def update_plot(self, frame):
        self.line.set_data(self.t[:frame], self.cooling_function(self.t[:frame]))
        self.ax.set_xlim(0, self.t[frame] if self.t[frame] > 10 else 10)
        return self.line,

    def reset_plot(self):
        self.anim.event_source.stop()
        self.line.set_data([], [])
        self.anim = FuncAnimation(self.fig, self.update_plot, frames=len(self.t),
                                  interval=30, blit=False, repeat=True)
        self.canvas.draw()

if __name__ == "__main__":
    root = tk.Tk()
    app = CoolingApp(root)
    root.mainloop()
