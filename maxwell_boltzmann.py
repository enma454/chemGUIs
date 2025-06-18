"""
@author: enmanuelduarte

GUI to understand the Maxwell-Boltzmann distribution of molecular speeds. 
The temperature effect on the speeds can be found along with its effect on lighter or heavier molecules. 

"""


import tkinter as tk
from tkinter import ttk, messagebox
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import scipy.constants as const
import logging

# Set up logging to console for diagnostics
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

# Maxwell-Boltzmann speed distribution function (M in g/mol)
def maxwell_distribution(v, M, T):
    R = const.R  # Gas constant (J/(mol·K), 8.314)
    try:
        # Prefactor: (M / (2 * pi * R * T))^(3/2)
        factor = (M / (2 * np.pi * R * T)) ** (3/2)
        if not np.isfinite(factor):
            logging.error("Prefactor is non-finite: M=%e g/mol, T=%e K, factor=%e", M, T, factor)
            return None
        # Exponential term: exp(-M * v^2 / (2 * R * T))
        exponent = -M * v**2 / (2 * R * T)
        if np.any(exponent < -700):  # Prevent underflow
            logging.warning("Exponent too negative, clipping: min exponent=%e", np.min(exponent))
            exponent = np.clip(exponent, -700, 0)
        prob = 4 * np.pi * v**2 * factor * np.exp(exponent)
        if not np.all(np.isfinite(prob)):
            logging.error("Probability contains non-finite values: min=%e, max=%e", np.min(prob), np.max(prob))
            return None
        return prob
    except Exception as e:
        logging.error("Error in distribution calculation: %s", str(e))
        return None

# Main GUI class
class MaxwellGUI:
    def __init__(self, root):
        self.root = root
        self.root.title("Maxwell-Boltzmann Speed Distribution")
        
        # Default values
        self.M = tk.DoubleVar(value=2.0)        # Molar mass in g/mol (H2)
        self.T = tk.DoubleVar(value=100.0)      # Temperature in K
        self.v_max = tk.DoubleVar(value=200.0)  # Max speed in m/s
        
        # Create input frame
        input_frame = ttk.LabelFrame(root, text="Parameters", padding=10)
        input_frame.grid(row=0, column=0, padx=10, pady=10, sticky="nsew")
        
        # Molar mass slider (2 to 100 g/mol)
        ttk.Label(input_frame, text="Molar Mass (g/mol):").grid(row=0, column=0, padx=5, pady=5)
        m_slider = ttk.Scale(input_frame, from_=2.0, to=100.0, orient="horizontal", 
                            variable=self.M, command=self.update_plot)
        m_slider.grid(row=0, column=1, padx=5, pady=5)
        self.m_label = ttk.Label(input_frame, text=f"{self.M.get():.1f}")
        self.m_label.grid(row=0, column=2, padx=5, pady=5)
        
        # Temperature slider (100 to 1000 K)
        ttk.Label(input_frame, text="Temperature (K):").grid(row=1, column=0, padx=5, pady=5)
        T_slider = ttk.Scale(input_frame, from_=100.0, to=1000.0, orient="horizontal", 
                            variable=self.T, command=self.update_plot)
        T_slider.grid(row=1, column=1, padx=5, pady=5)
        self.T_label = ttk.Label(input_frame, text=f"{self.T.get():.0f}")
        self.T_label.grid(row=1, column=2, padx=5, pady=5)
        
        # Speed range input
        ttk.Label(input_frame, text="Max Speed (m/s):").grid(row=2, column=0, padx=5, pady=5)
        ttk.Entry(input_frame, textvariable=self.v_max).grid(row=2, column=1, padx=5, pady=5)
        
        # Update button for speed range
        ttk.Button(input_frame, text="Update Plot", command=self.update_plot).grid(row=3, column=0, columnspan=2, pady=10)
        
        # Matplotlib plot
        self.fig, self.ax = plt.subplots(figsize=(6, 4))
        self.canvas = FigureCanvasTkAgg(self.fig, master=root)
        self.canvas.get_tk_widget().grid(row=0, column=1, padx=10, pady=10)
        
        # Initial plot
        self.update_plot()
    
    def update_plot(self, *args):
        # Update slider labels
        self.m_label.config(text=f"{self.M.get():.1f}")
        self.T_label.config(text=f"{self.T.get():.0f}")
        
        # Get parameters with validation
        try:
            M = self.M.get()  # Molar mass in g/mol
            T = self.T.get()
            v_max = self.v_max.get()
            if M <= 0 or T <= 0 or v_max <= 0:
                messagebox.showerror("Invalid Input", "Parameters must be positive.")
                return
            if v_max > 5000:
                messagebox.showerror("Invalid Input", "Max speed should be <= 5000 m/s.")
                return
        except tk.TclError:
            messagebox.showerror("Invalid Input", "Please enter valid numbers.")
            return
        
        # Speed points
        v = np.linspace(1e-6, v_max, 1000)
        
        # Calculate distribution
        prob = maxwell_distribution(v, M, T)
        if prob is None:
            messagebox.showerror("Calculation Error", "Distribution calculation failed. Try different parameters (e.g., lower molar mass or higher temperature).")
            return
        if np.max(prob) == 0:
            messagebox.showerror("Calculation Error", "Distribution is zero everywhere. Try increasing temperature or decreasing molar mass.")
            return
        
        # Plot
        self.ax.clear()
        self.ax.plot(v, prob, label='Maxwell-Boltzmann', color='blue')
        self.ax.set_xlabel('Speed (m/s)')
        self.ax.set_ylabel('Probability Density (a.u.)')
        self.ax.set_title('Maxwell-Boltzmann Speed Distribution')
        self.ax.legend()
        self.ax.grid(True)
        self.canvas.draw()
        
        # Log peak speed for verification
        peak_idx = np.argmax(prob)
        peak_v = v[peak_idx]
        logging.info("Plot updated: M=%f g/mol, T=%f K, v_max=%f m/s, peak speed=%f m/s", M, T, v_max, peak_v)

# Run the application
if __name__ == "__main__":
    root = tk.Tk()
    app = MaxwellGUI(root)
    root.mainloop()