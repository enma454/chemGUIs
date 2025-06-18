"""
@author: enmanuelduarte

GUI to explore the chemical kinetics for the mechanism A -> B -> C. 
By manipulating the parameters (initial concentrations and kinetic constants) 
it is possible to explore distinct types of behavior of chemical species in relation to reactivity.
"""

import tkinter as tk
from tkinter import ttk
import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg

# Define the differential equations for A -> B -> C
def kinetic_model(y, t, k1, k2):
    A, B, C = y
    dA_dt = -k1 * A
    dB_dt = k1 * A - k2 * B
    dC_dt = k2 * B
    return [dA_dt, dB_dt, dC_dt]

# Main GUI class
class KineticGUI:
    def __init__(self, root):
        self.root = root
        self.root.title("Kinetic Mechanism A -> B -> C")
        
        # Default values
        self.A0 = tk.DoubleVar(value=1.0)
        self.B0 = tk.DoubleVar(value=0.0)
        self.C0 = tk.DoubleVar(value=0.0)
        self.t_max = tk.DoubleVar(value=10.0)
        self.k1 = tk.DoubleVar(value=1.0)
        self.k2 = tk.DoubleVar(value=1.0)
        
        # Create input frame
        input_frame = ttk.LabelFrame(root, text="Parameters", padding=10)
        input_frame.grid(row=0, column=0, padx=10, pady=10, sticky="nsew")
        
        # Concentration inputs
        ttk.Label(input_frame, text="Initial A (M):").grid(row=0, column=0, padx=5, pady=5)
        ttk.Entry(input_frame, textvariable=self.A0).grid(row=0, column=1, padx=5, pady=5)
        
        ttk.Label(input_frame, text="Initial B (M):").grid(row=1, column=0, padx=5, pady=5)
        ttk.Entry(input_frame, textvariable=self.B0).grid(row=1, column=1, padx=5, pady=5)
        
        ttk.Label(input_frame, text="Initial C (M):").grid(row=2, column=0, padx=5, pady=5)
        ttk.Entry(input_frame, textvariable=self.C0).grid(row=2, column=1, padx=5, pady=5)
        
        ttk.Label(input_frame, text="Time Range (s):").grid(row=3, column=0, padx=5, pady=5)
        ttk.Entry(input_frame, textvariable=self.t_max).grid(row=3, column=1, padx=5, pady=5)
        
        # Kinetic constant sliders
        ttk.Label(input_frame, text="k1 (1/s):").grid(row=4, column=0, padx=5, pady=5)
        k1_slider = ttk.Scale(input_frame, from_=0.1, to=5.0, orient="horizontal", 
                            variable=self.k1, command=self.update_plot)
        k1_slider.grid(row=4, column=1, padx=5, pady=5)
        self.k1_label = ttk.Label(input_frame, text=f"{self.k1.get():.2f}")
        self.k1_label.grid(row=4, column=2, padx=5, pady=5)
        
        ttk.Label(input_frame, text="k2 (1/s):").grid(row=5, column=0, padx=5, pady=5)
        k2_slider = ttk.Scale(input_frame, from_=0.1, to=5.0, orient="horizontal", 
                            variable=self.k2, command=self.update_plot)
        k2_slider.grid(row=5, column=1, padx=5, pady=5)
        self.k2_label = ttk.Label(input_frame, text=f"{self.k2.get():.2f}")
        self.k2_label.grid(row=5, column=2, padx=5, pady=5)
        
        # Update button for concentration and time inputs
        ttk.Button(input_frame, text="Update Plot", command=self.update_plot).grid(row=6, column=0, columnspan=2, pady=10)
        
        # Matplotlib plot
        self.fig, self.ax = plt.subplots(figsize=(6, 4))
        self.canvas = FigureCanvasTkAgg(self.fig, master=root)
        self.canvas.get_tk_widget().grid(row=0, column=1, padx=10, pady=10)
        
        # Initial plot
        self.update_plot()
    
    def update_plot(self, *args):
        # Update slider labels
        self.k1_label.config(text=f"{self.k1.get():.2f}")
        self.k2_label.config(text=f"{self.k2.get():.2f}")
        
        # Get parameters
        try:
            A0 = self.A0.get()
            B0 = self.B0.get()
            C0 = self.C0.get()
            t_max = self.t_max.get()
            k1 = self.k1.get()
            k2 = self.k2.get()
        except tk.TclError:
            return  # Invalid input, skip update
        
        # Time points
        t = np.linspace(0, t_max, 1000)
        
        # Solve ODE
        y0 = [A0, B0, C0]
        solution = odeint(kinetic_model, y0, t, args=(k1, k2))
        
        # Plot
        self.ax.clear()
        self.ax.plot(t, solution[:, 0], label='A', color='blue')
        self.ax.plot(t, solution[:, 1], label='B', color='green')
        self.ax.plot(t, solution[:, 2], label='C', color='red')
        self.ax.set_xlabel('Time (s)')
        self.ax.set_ylabel('Concentration (M)')
        self.ax.set_title('A -> B -> C Kinetic Mechanism')
        self.ax.legend()
        self.ax.grid(True)
        self.canvas.draw()

# Run the application
if __name__ == "__main__":
    root = tk.Tk()
    app = KineticGUI(root)
    root.mainloop()