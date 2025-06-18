"""
@author: enmanuelduarte

GUI for understanding the relationship between kinetics and thermodynamics (chemical equilibrium) in chemical reactions. 
This relationship can be explored by changing the kinetic constants k1, k2, the initial concentrations of reactant A0 and product B0, and the temperature T.

"""

import tkinter as tk
from tkinter import ttk, messagebox
import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk

class ChemicalEquilibriumApp:
    def __init__(self, root):
        self.root = root
        self.root.title("Chemical Equilibrium Simulator: A <=> B")

        # Constants
        self.R = 8.314  # Gas constant in J/(mol·K)

        # --- GUI Setup ---
        self.create_parameter_frame()
        self.create_plot_frame()
        self.update_plots() # Initial plot

    def create_parameter_frame(self):
        """Creates the frame for input parameters (sliders and text entries)."""
        param_frame = ttk.LabelFrame(self.root, text="Parameters")
        param_frame.grid(row=0, column=0, padx=10, pady=10, sticky="nsew")

        # K1 (rate constant)
        ttk.Label(param_frame, text="k1 (0.01 - 1):").grid(row=0, column=0, padx=5, pady=5, sticky="w")
        self.k1_var = tk.DoubleVar(value=0.4)
        self.k1_slider = ttk.Scale(param_frame, from_=0.01, to=1.0, orient="horizontal",
                                   variable=self.k1_var, command=self.update_plots_from_slider)
        self.k1_slider.grid(row=0, column=1, padx=5, pady=5, sticky="ew")
        self.k1_label = ttk.Label(param_frame, text=f"{self.k1_var.get():.2f}")
        self.k1_label.grid(row=0, column=2, padx=5, pady=5, sticky="w")

        # K2 (rate constant)
        ttk.Label(param_frame, text="k2 (0.01 - 1):").grid(row=1, column=0, padx=5, pady=5, sticky="w")
        self.k2_var = tk.DoubleVar(value=0.2)
        self.k2_slider = ttk.Scale(param_frame, from_=0.01, to=1.0, orient="horizontal",
                                   variable=self.k2_var, command=self.update_plots_from_slider)
        self.k2_slider.grid(row=1, column=1, padx=5, pady=5, sticky="ew")
        self.k2_label = ttk.Label(param_frame, text=f"{self.k2_var.get():.2f}")
        self.k2_label.grid(row=1, column=2, padx=5, pady=5, sticky="w")

        # A0 (initial concentration)
        ttk.Label(param_frame, text="A0 (0 - 10):").grid(row=2, column=0, padx=5, pady=5, sticky="w")
        self.A0_entry = ttk.Entry(param_frame, width=10)
        self.A0_entry.insert(0, "1.0")
        self.A0_entry.grid(row=2, column=1, padx=5, pady=5, sticky="w")

        # B0 (initial concentration)
        ttk.Label(param_frame, text="B0 (0 - 10):").grid(row=3, column=0, padx=5, pady=5, sticky="w")
        self.B0_entry = ttk.Entry(param_frame, width=10)
        self.B0_entry.insert(0, "0.0")
        self.B0_entry.grid(row=3, column=1, padx=5, pady=5, sticky="w")

        # Temperature (T)
        ttk.Label(param_frame, text="T (200 - 500 K):").grid(row=4, column=0, padx=5, pady=5, sticky="w")
        self.T_var = tk.DoubleVar(value=298.15)
        self.T_slider = ttk.Scale(param_frame, from_=200, to=500, orient="horizontal",
                                  variable=self.T_var, command=self.update_plots_from_slider)
        self.T_slider.grid(row=4, column=1, padx=5, pady=5, sticky="ew")
        self.T_label = ttk.Label(param_frame, text=f"{self.T_var.get():.1f} K")
        self.T_label.grid(row=4, column=2, padx=5, pady=5, sticky="w")

        # Plot Button
        self.plot_button = ttk.Button(param_frame, text="Update Plots", command=self.update_plots)
        self.plot_button.grid(row=5, column=0, columnspan=3, pady=10)

        # Information Display
        self.info_label = ttk.Label(param_frame, text="Initial Q / K info will appear here.", wraplength=250, justify="left")
        self.info_label.grid(row=6, column=0, columnspan=3, padx=5, pady=5, sticky="w")

        # Configure column weights to make sliders expand
        param_frame.grid_columnconfigure(1, weight=1)


    def create_plot_frame(self):
        """Creates the frame to hold the matplotlib plots."""
        plot_frame = ttk.Frame(self.root)
        plot_frame.grid(row=0, column=1, padx=10, pady=10, sticky="nsew")
        self.root.grid_columnconfigure(1, weight=1)
        self.root.grid_rowconfigure(0, weight=1)

        self.fig, (self.ax1, self.ax2) = plt.subplots(2, 1, figsize=(8, 7), layout='constrained')
        self.canvas = FigureCanvasTkAgg(self.fig, master=plot_frame)
        self.canvas_widget = self.canvas.get_tk_widget()
        self.canvas_widget.pack(side=tk.TOP, fill=tk.BOTH, expand=True)

        self.toolbar = NavigationToolbar2Tk(self.canvas, plot_frame)
        self.toolbar.update()
        self.canvas_widget.pack(side=tk.TOP, fill=tk.BOTH, expand=True)

    def update_plots_from_slider(self, event):
        """Updates the numerical labels next to sliders and triggers plot update."""
        self.k1_label.config(text=f"{self.k1_var.get():.2f}")
        self.k2_label.config(text=f"{self.k2_var.get():.2f}")
        self.T_label.config(text=f"{self.T_var.get():.1f} K")
        self.update_plots()

    def update_plots(self):
        """Fetches parameters, performs calculations, and updates both plots."""
        try:
            k1 = self.k1_var.get()
            k2 = self.k2_var.get()
            A0 = float(self.A0_entry.get())
            B0 = float(self.B0_entry.get())
            T = self.T_var.get()

            # Input validation for A0, B0
            if not (0 <= A0 <= 10 and 0 <= B0 <= 10):
                messagebox.showerror("Input Error", "Initial concentrations A0 and B0 must be between 0 and 10.")
                return
            if A0 <= 0 and B0 <= 0:
                 messagebox.showerror("Input Error", "At least one initial concentration (A0 or B0) must be greater than 0.")
                 return

            # Calculate Equilibrium Constant K
            if k2 == 0: # Avoid division by zero
                messagebox.showwarning("Warning", "k2 cannot be zero. Setting it to a small value for calculation.")
                k2 = 0.0001 # Set to a very small non-zero value
                self.k2_var.set(k2) # Update slider
                self.k2_label.config(text=f"{self.k2_var.get():.4f}") # Update label

            K_eq = k1 / k2

            # Calculate Initial Reaction Quotient Q
            # If A0 is 0, Q is infinity (if B0 > 0), if B0 is 0, Q is 0 (if A0 > 0)
            if A0 <= 1e-9 and B0 > 0: # A0 is practically zero, B0 is positive
                Q_initial = np.inf
            elif B0 <= 1e-9 and A0 > 0: # B0 is practically zero, A0 is positive
                Q_initial = 0.0
            elif A0 <= 1e-9 and B0 <= 1e-9: # Both are practically zero (should be caught by earlier validation)
                 messagebox.showerror("Error", "Both A0 and B0 are close to zero; cannot calculate Q.")
                 return
            else:
                Q_initial = B0 / A0

            # Calculate equilibrium extent of reaction (xi_eq)
            # Total initial concentration
            total_conc = A0 + B0
            if K_eq <= 0: # Handle cases where K_eq is non-positive
                A_eq_val = total_conc # Effectively all A
            else:
                A_eq_val = total_conc / (1 + K_eq) # Equilibrium [A]

            xi_eq = A0 - A_eq_val # Extent of reaction at equilibrium

            # Update Q/K information label
            info_text = f"Equilibrium Constant (K) = {K_eq:.2f}\n"
            info_text += f"Initial Reaction Quotient (Q) = {Q_initial:.2f}\n"
            if Q_initial > K_eq:
                info_text += "Q > K: Reaction favors reactants (shifts left)\n"
            elif Q_initial < K_eq:
                info_text += "Q < K: Reaction favors products (shifts right)\n"
            else:
                info_text += "Q = K: System is at equilibrium\n"

            info_text += f"Equilibrium Extent of Reaction (ξ_eq) = {xi_eq:.2f}\n"
            if not (0 <= xi_eq <= 1):
                info_text += "(Note: ξ_eq is outside the plotted range 0-1)"

            self.info_label.config(text=info_text)

            # --- Plot 1: Concentration vs. Time ---
            self.plot_concentration_vs_time(k1, k2, A0, B0)

            # --- Plot 2: Gibbs Free Energy vs. Extent of Reaction ---
            self.plot_gibbs_vs_extent(k1, k2, A0, B0, T, K_eq)

            self.canvas.draw()

        except ValueError:
            messagebox.showerror("Input Error", "Please enter valid numerical values for initial concentrations.")
        except Exception as e:
            messagebox.showerror("An Error Occurred", f"An unexpected error occurred: {e}")

    def reaction_odes(self, C, t, k1, k2):
        """
        Defines the system of ordinary differential equations for the reaction A <=> B.
        C[0] = [A], C[1] = [B]
        """
        A_conc, B_conc = C
        dAdt = -k1 * A_conc + k2 * B_conc
        dBdt = k1 * A_conc - k2 * B_conc
        return [dAdt, dBdt]

    def plot_concentration_vs_time(self, k1, k2, A0, B0):
        """Calculates and plots concentration vs. time."""
        t_span = np.linspace(0, 20, 500)  # Time from 0 to 20 units
        initial_concentrations = [A0, B0]

        # Solve ODEs
        sol = odeint(self.reaction_odes, initial_concentrations, t_span, args=(k1, k2))
        A_conc = sol[:, 0]
        B_conc = sol[:, 1]

        self.ax1.clear()
        self.ax1.plot(t_span, A_conc, label='[A]/M')
        self.ax1.plot(t_span, B_conc, label='[B]/M')
        self.ax1.set_xlabel('Time / s')
        self.ax1.set_ylabel('Concentration / M')
        self.ax1.set_title('Concentration vs. Time')
        self.ax1.legend()
        self.ax1.grid(True)
        self.ax1.set_ylim(bottom=0) # Ensure y-axis starts at 0

    def plot_gibbs_vs_extent(self, k1, k2, A0, B0, T, K_eq):
        """Calculates and plots Gibbs free energy of the system vs. extent of reaction."""
        # Calculate standard Gibbs free energy change
        if K_eq <= 1e-9: # If K_eq is very small, delta_G_std_0 is very positive (favors reactants)
            delta_G_std_0 = 100000.0 # Assign a large positive value
        elif K_eq >= 1e9: # If K_eq is very large, delta_G_std_0 is very negative (favors products)
            delta_G_std_0 = -100000.0 # Assign a large negative value
        else:
            delta_G_std_0 = -self.R * T * np.log(K_eq)

        # Extent of reaction (xi) always from 0 to 1
        xi_values = np.linspace(0, 1.0, 100)

        gibbs_system_values = []
        for xi in xi_values:
            # Ensure concentrations remain positive to avoid log errors
            A_at_xi = max(1e-9, A0 - xi)
            B_at_xi = max(1e-9, B0 + xi)

            # Calculate total Gibbs free energy of the system
            # G_sys = (A0 - xi)*mu_A + (B0 + xi)*mu_B
            # where mu_i = mu_i_standard + RT ln C_i
            # We can use a relative G_sys to show the parabolic shape:
            # G_relative = (B0 + xi)*delta_G_standard + RT*((A0 - xi)ln(A0 - xi) + (B0 + xi)ln(B0 + xi))
            
            # The calculation used here is based on the total Gibbs free energy relative to a reference.
            # G_system_plot = -RT * (B0 + xi) * ln(K_eq) + RT * ( (A0 - xi) * ln(A0 - xi) + (B0 + xi) * ln(B0 + xi) )
            # This is equivalent to RT * [ (A_0 - xi)ln(A_0 - xi) + (B_0 + xi)ln(B_0 + xi) + (xi)ln(K_eq) - (A_0)ln(A_0) - (B_0)ln(B_0)]
            # which is (G_system - G_system_initial) - constant term related to standard states.
            # Use the explicit form for G_plot(xi) derived for the U-shape.
            
            term1 = -(B0 + xi) * np.log(K_eq)
            term2 = (A0 - xi) * np.log(A0 - xi)
            term3 = (B0 + xi) * np.log(B0 + xi)

            G_plot = self.R * T * (term1 + term2 + term3)
            gibbs_system_values.append(G_plot)


        self.ax2.clear()
        self.ax2.plot(xi_values, gibbs_system_values, label='G_system vs. Extent of Reaction')
        self.ax2.set_xlabel('Extent of Reaction (ξ)')
        self.ax2.set_ylabel('Gibbs Free Energy (J/mol)')
        self.ax2.set_title('Gibbs Free Energy of the System vs. Extent of Reaction')
        # The minimum of this curve represents equilibrium; ΔG=0 is the slope at this minimum.
        self.ax2.legend()
        self.ax2.grid(True)
        # Set y-limits dynamically to show the relevant part of the curve
        if gibbs_system_values:
            finite_gibbs_values = [val for val in gibbs_system_values if np.isfinite(val)]
            if finite_gibbs_values:
                min_g = min(finite_gibbs_values)
                max_g = max(finite_gibbs_values)
                padding = (max_g - min_g) * 0.1
                self.ax2.set_ylim(min_g - padding, max_g + padding)
            else: # Fallback if all values are infinite (e.g., extremely high/low K_eq)
                self.ax2.set_ylim(-20000, 20000)
        else:
            self.ax2.set_ylim(-20000, 20000)


# Main application entry point
if __name__ == "__main__":
    root = tk.Tk()
    app = ChemicalEquilibriumApp(root)
    root.mainloop()