import numpy as np
from Solver.Processor import Processor
from Solver.Post_Processor import Post_Processor as post
from Solver.Optimization import Optimizer
from Config.config import load_config

import os

dir = "Results"

class Run:
    def __init__(self):
        self.config = load_config()
        self._load_parameters()

    def _load_parameters(self):
        """Extract parameters from the config file."""
        impedance_config = self.config["impedance"]

        self.frequencies = np.linspace(
            impedance_config["frequencies"]["start"],
            impedance_config["frequencies"]["end"],
            impedance_config["frequencies"]["num_points"]
        )
        self.altitude = impedance_config["altitude"]
        
        self.p_acous_pa = self.config["p_acous_pa"]

        # Liner's geometry
        liner = self.config["liner"]
        self.L, self.d, self.sigma, self.e = liner["L"], liner["d"], liner["sigma"], liner["e"]

    def run_impedance_study(self):
        """Run impedance calculations and plot results."""
        varying_params = {
            "L": np.linspace(10e-3, 20e-3, 5),
            "d": np.linspace(1e-3, 2.e-3, 5),
            "sigma": np.linspace(0.1, 0.2, 5),
            "e": np.linspace(1e-3, 2.e-3, 5),
        }

        for param, values in varying_params.items():
            if param == "L":
                impedances = Processor.compute_impedance_varying_param(
                    param, values, self.d, self.sigma, self.e, self.altitude, self.frequencies, self.p_acous_pa
                )
            elif param == "d":
                impedances = Processor.compute_impedance_varying_param(
                    param, self.L, values, self.sigma, self.e, self.altitude, self.frequencies, self.p_acous_pa
                )
            elif param == "sigma":
                impedances = Processor.compute_impedance_varying_param(
                    param, self.L, self.d, values, self.e, self.altitude, self.frequencies, self.p_acous_pa
                )
            elif param == "e":
                impedances = Processor.compute_impedance_varying_param(
                    param, self.L, self.d, self.sigma, values, self.altitude, self.frequencies, self.p_acous_pa
                )
            save_path = f"{dir}/impedance_study/{param}_absorption.png"
            post.plot_absorption_coefficients(self.frequencies, impedances, save_path)

    def run_optimization(self):
        """Run optimization to find the best liner geometry."""
        optimization = self.config["optimization"]
        frequencie = optimization["frequencie"]
        altitudes = np.linspace(
            optimization["altitudes"]["start"],
            optimization["altitudes"]["end"],
            optimization["altitudes"]["num_points"]
        )

        L_var = np.linspace(*optimization["varying_params"]["L_range"])
        sigma_var = np.linspace(*optimization["varying_params"]["sigma_range"])

        alpha = Optimizer.loss_function_values(L_var, self.d, sigma_var, self.e, self.altitude, frequencie, self.p_acous_pa)
        alpha_max, optimum_L, optimum_sigma = Optimizer.optimum_geometry(L_var, self.d, sigma_var, self.e, altitudes, frequencie, self.p_acous_pa)
        
        save_path = f"{dir}/Geometry_study/"

        post.plot_loss_function(sigma_var, L_var, alpha, save_path+"Alpha_mapping.png")
        post.plot_optimum_geometry(altitudes, alpha_max, optimum_L, optimum_sigma, save_path+"Optimum_geometry.png")


    def run_optimization2(self):
        """Run optimization to find the best liner geometry."""
        optimization = self.config["optimization"]
        frequencie = optimization["frequencie"]
        altitudes = np.linspace(
            optimization["altitudes"]["start"],
            optimization["altitudes"]["end"],
            optimization["altitudes"]["num_points"]
        )

        L_var = np.linspace(*optimization["varying_params"]["L_range"])
        sigma_var = np.linspace(*optimization["varying_params"]["sigma_range"])

        for i,altitude in enumerate(altitudes):
            alpha = Optimizer.loss_function_values(L_var, self.d, sigma_var, self.e, altitude, frequencie, self.p_acous_pa)
       
            save_path = f"{dir}/Geometry_study/"

            post.plot_loss_function(sigma_var, L_var, alpha, save_path+ str(i)+".png")
