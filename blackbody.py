"""
An interactive GUI to explore the behavior of the blackbody radiation spectrum on temperature changes. 

@author: enmanuelduarte

"""

import tkinter as tk
from tkinter import ttk
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import numpy as np

class BlackbodyRadiationExplorer:
    def __init__(self, root):
        self.root = root
        self.root.title("Blackbody Radiation Explorer")
        self.root.geometry("900x700")
        
        # Physical constants
        self.h = 6.626e-34  # Planck's constant (J⋅s)
        self.c = 3e8        # Speed of light (m/s)
        self.k = 1.381e-23  # Boltzmann constant (J/K)
        
        # Temperature range and initial value
        self.temp_min = 500
        self.temp_max = 8000
        self.initial_temp = 3000
        
        self.setup_gui()
        self.update_plot()
        
    def setup_gui(self):
        # Create main frame
        main_frame = ttk.Frame(self.root, padding="10")
        main_frame.grid(row=0, column=0, sticky=(tk.W, tk.E, tk.N, tk.S))
        
        # Configure grid weights
        self.root.columnconfigure(0, weight=1)
        self.root.rowconfigure(0, weight=1)
        main_frame.columnconfigure(0, weight=1)
        main_frame.rowconfigure(1, weight=1)
        
        # Control frame
        control_frame = ttk.Frame(main_frame)
        control_frame.grid(row=0, column=0, sticky=(tk.W, tk.E), pady=(0, 10))
        control_frame.columnconfigure(1, weight=1)
        
        # Temperature slider
        ttk.Label(control_frame, text="Temperature (K):").grid(row=0, column=0, padx=(0, 10))
        
        self.temp_var = tk.DoubleVar(value=self.initial_temp)
        self.temp_slider = ttk.Scale(
            control_frame, 
            from_=self.temp_min, 
            to=self.temp_max,
            variable=self.temp_var,
            orient=tk.HORIZONTAL,
            command=self.on_temperature_change
        )
        self.temp_slider.grid(row=0, column=1, sticky=(tk.W, tk.E), padx=(0, 10))
        
        # Temperature value label
        self.temp_label = ttk.Label(control_frame, text=f"{self.initial_temp:.0f} K")
        self.temp_label.grid(row=0, column=2)
        
        # Information frame
        info_frame = ttk.Frame(main_frame)
        info_frame.grid(row=0, column=1, sticky=(tk.W, tk.E), padx=(20, 0))
        
        # Peak wavelength info
        self.peak_info = ttk.Label(info_frame, text="", font=("Arial", 10))
        self.peak_info.grid(row=0, column=0, sticky=tk.W)
        
        # Wien's law info
        self.wien_info = ttk.Label(info_frame, text="Wien's Law: λ_max = 2898 μm⋅K / T", 
                                  font=("Arial", 9), foreground="gray")
        self.wien_info.grid(row=1, column=0, sticky=tk.W)
        
        # Create matplotlib figure
        self.fig, self.ax = plt.subplots(figsize=(10, 6))
        self.fig.tight_layout(pad=3.0)
        
        # Embed plot in tkinter
        plot_frame = ttk.Frame(main_frame)
        plot_frame.grid(row=1, column=0, columnspan=2, sticky=(tk.W, tk.E, tk.N, tk.S))
        plot_frame.columnconfigure(0, weight=1)
        plot_frame.rowconfigure(0, weight=1)
        
        self.canvas = FigureCanvasTkAgg(self.fig, plot_frame)
        self.canvas.get_tk_widget().grid(row=0, column=0, sticky=(tk.W, tk.E, tk.N, tk.S))
        
        # Add preset temperature buttons
        preset_frame = ttk.Frame(main_frame)
        preset_frame.grid(row=2, column=0, columnspan=2, pady=(10, 0))
        
        ttk.Label(preset_frame, text="Preset Temperatures:").grid(row=0, column=0, padx=(0, 10))
        
        presets = [
            ("Sun Surface (~5778 K)", 5778),
            ("Incandescent Bulb (~2700 K)", 2700),
            ("Human Body (~310 K)", 310),
            ("Room Temp (~300 K)", 300)
        ]
        
        for i, (label, temp) in enumerate(presets):
            btn = ttk.Button(preset_frame, text=label, 
                           command=lambda t=temp: self.set_temperature(t))
            btn.grid(row=0, column=i+1, padx=5)
    
    def planck_function(self, wavelength_m, temperature):
        """
        Calculate Planck's blackbody radiation formula
        wavelength_m: wavelength in meters
        temperature: temperature in Kelvin
        Returns: spectral radiance in W⋅sr⁻¹⋅m⁻³
        """
        # Avoid division by zero and overflow
        with np.errstate(over='ignore', divide='ignore', invalid='ignore'):
            numerator = 2 * self.h * self.c**2 / wavelength_m**5
            denominator = np.exp(self.h * self.c / (wavelength_m * self.k * temperature)) - 1
            
            # Handle potential overflow/underflow
            result = np.where(denominator > 0, numerator / denominator, 0)
            
        return result
    
    def update_plot(self):
        temperature = self.temp_var.get()
        
        # Fixed wavelength range to show all curves on same scale
        wavelength_microns = np.linspace(0.1, 20, 1000)
        wavelength_m = wavelength_microns * 1e-6  # Convert to meters
        
        # Calculate Planck function (actual energy density)
        spectral_radiance = self.planck_function(wavelength_m, temperature)
        
        # Scale to make values more manageable for display (but preserve relative scaling)
        # Use a fixed reference temperature to maintain relative scaling between curves
        reference_temp = 1000  # K
        reference_peak = self.planck_function(np.array([2898e-6 / reference_temp]), reference_temp)[0]
        
        # Scale all curves relative to this reference
        scaled_energy = spectral_radiance / reference_peak
        
        # Clear and plot
        self.ax.clear()
        self.ax.plot(wavelength_microns, scaled_energy, 'b-', linewidth=2, 
                    label=f'T = {temperature:.0f} K')
        
        # Find and mark the peak
        wien_peak = 2898 / temperature  # μm
        if wien_peak >= wavelength_microns.min() and wien_peak <= wavelength_microns.max():
            peak_idx = np.argmin(np.abs(wavelength_microns - wien_peak))
            peak_energy = scaled_energy[peak_idx]
            self.ax.plot(wien_peak, peak_energy, 'ro', markersize=8, 
                        label=f'Peak: {wien_peak:.2f} μm')
        
        # Styling
        self.ax.set_xlabel('Wavelength (μm)', fontsize=12)
        self.ax.set_ylabel('Energy Density (relative to 1000K reference)', fontsize=12)
        self.ax.set_title(f'Blackbody Radiation Spectrum (T = {temperature:.0f} K)', fontsize=14)
        self.ax.grid(True, alpha=0.3)
        
        # Fixed axis limits to show the physics clearly
        self.ax.set_xlim(0.1, 20)
        
        # Dynamic y-limit based on temperature to show the height changes
        max_energy = np.max(scaled_energy)
        self.ax.set_ylim(0, max_energy * 1.1)
        
        # Add visible light region 
        self.ax.axvspan(0.38, 0.75, alpha=0.2, color='yellow', label='Visible Light')
        
        # Add UV and infrared regions
        self.ax.axvspan(0.1, 0.38, alpha=0.1, color='purple', label='UV')
        self.ax.axvspan(0.75, 20, alpha=0.1, color='red', label='Infrared')
        
        self.ax.legend(loc='upper right')
        
        # Update canvas
        self.canvas.draw()
        
        # Update info labels
        self.update_info_labels(temperature, wien_peak, max_energy)
    
    def update_info_labels(self, temperature, wien_peak, max_energy):
        self.temp_label.config(text=f"{temperature:.0f} K")
        
        # Peak wavelength info
        peak_text = f"Peak wavelength: {wien_peak:.2f} μm"
        
        # Determine which region the peak is in
        if wien_peak < 0.38:
            region = "UV"
        elif wien_peak <= 0.75:
            region = "Visible"
        else:
            region = "IR"
        
        peak_text += f" ({region})"
        
        # Add energy information
        energy_text = f"\nPeak energy density: {max_energy:.2f}x reference"
        
        # Stefan-Boltzmann total power scaling (T^4)
        reference_temp = 1000
        total_power_ratio = (temperature / reference_temp) ** 4
        energy_text += f"\nTotal radiated power: {total_power_ratio:.2f}x reference"
        
        self.peak_info.config(text=peak_text + energy_text)
    
    def on_temperature_change(self, value):
        self.update_plot()
    
    def set_temperature(self, temp):
        # Clamp temperature to slider range
        temp = max(self.temp_min, min(self.temp_max, temp))
        self.temp_var.set(temp)
        self.update_plot()

def main():
    root = tk.Tk()
    app = BlackbodyRadiationExplorer(root)
    root.mainloop()

if __name__ == "__main__":
    main()