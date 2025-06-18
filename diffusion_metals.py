"""
An interactive GUI to shed light on a diffusion couple formed by joining two pure metal bars together 
that present an intimate contact between the two metal surfaces. 
The couple is heated for an extended period of time at an elevated temperature 
(but below the melting temperatures of both metals) and cooled to room temperature. 
The concentration of metals vs position in time can be found by modifying the diffusion coefficients of both metals. 

@author: enmanuelduarte

"""
import tkinter as tk
from tkinter import ttk
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from scipy.special import erf

class DiffusionGUI:
    def __init__(self, root):
        self.root = root
        self.root.title("Diffusion Couple Simulator")

        # Default values
        self.D1 = tk.DoubleVar(value=1e-14) # Diffusivity of metal A (m2/s)
        self.D2 = tk.DoubleVar(value=1e-14) # Diffusivity of metal B (m2/s)

        # Layout
        self.create_widgets()
        self.plot_graphs()

    def create_widgets(self):
        frame = ttk.Frame(self.root)
        frame.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)

        ttk.Label(frame, text="D1 (m²/s)").pack()
        ttk.Scale(frame, from_=1e-16, to=1e-10, variable=self.D1, orient='horizontal', command=self.on_change).pack(fill='x')

        ttk.Label(frame, text="D2 (m²/s)").pack()
        ttk.Scale(frame, from_=1e-16, to=1e-10, variable=self.D2, orient='horizontal', command=self.on_change).pack(fill='x')

        # Matplotlib figure
        self.fig, self.ax1 = plt.subplots(figsize=(6, 4))
        self.canvas = FigureCanvasTkAgg(self.fig, master=self.root)
        self.canvas.get_tk_widget().pack(side=tk.RIGHT, fill=tk.BOTH, expand=True)

    def concentration_profile(self, D, t):
        x = self.x
        return 0.5 * (1 + erf(x / (2 * np.sqrt(D * t))))

    def plot_graphs(self):
        self.x = np.linspace(-1e-3, 1e-3, 500)
        self.times = [10, 100, 1000, 5000, 10000]

        D1 = self.D1.get()
        D2 = self.D2.get()
        Deff = (D1 + D2) / 2

        self.ax1.clear()

        for t in self.times:
            conc = self.concentration_profile(Deff, t)
            self.ax1.plot(self.x * 1e6, conc, label=f't / s= {t:.0f}')

        self.ax1.set_title('Concentration Profile vs Position')
        self.ax1.set_xlabel('Position (µm)')
        self.ax1.set_ylabel('Concentration of Metal A / M')
        self.ax1.legend()
        self.ax1.grid(True)

        self.canvas.draw()

    def on_change(self, event):
        self.plot_graphs()

if __name__ == "__main__":
    root = tk.Tk()
    app = DiffusionGUI(root)
    root.mainloop()
