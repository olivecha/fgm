import numpy as np


class PickleFreeFlame(object):
    """
    A class to converting cantera FreeFlame
    objects to picklable objects
    """
    air =  {"O2":0.21, "N2":0.79}
    
    def __init__(self, source, phi, fuel, air=None):
        """
        Class constructor storing FreeFlame attributes
        as attributes of the current class
        """
        # Argument properties
        self.phi = phi
        self.fuel = fuel
        if air is not None:
            self.air = air
        
        # Base properties
        self.gas = source.gas
        self.grid = source.grid
        self.T = source.T
        self.P = source.P
        self.Y = source.Y
        self.X = source.X
        self.velocity = source.velocity

        # Elemental fractions
        self.elemental_mass_fraction = {}
        for e in ['H', 'C', 'N']:
            self.elemental_mass_fraction[e] = source.elemental_mass_fraction(e)
        self.elemental_mole_fraction = {}
        for e in ['H', 'C', 'N']:
            self.elemental_mole_fraction[e] = source.elemental_mole_fraction(e)

        # Physical properties
        self.density_mass = source.density_mass
        self.density = source.density_mass
        self.density_mole = source.density_mole
        self.volume_mass = source.volume_mass
        self.volume_mole = source.volume_mole
        self.int_energy_mass = source.int_energy_mass
        self.int_energy_mole = source.int_energy_mole
        self.enthalpy_mass = source.enthalpy_mass
        self.enthalpy_mole = source.enthalpy_mole
        self.entropy_mass = source.entropy_mass
        self.entropy_mole = source.entropy_mole
        self.cv_mass = source.cv_mass
        self.cv_mole = source.cv_mole
        self.cp_mass = source.cp_mass
        self.cp_mole = source.cp_mole
        self.viscosity = source.viscosity
        self.thermal_conductivity = source.thermal_conductivity
        self.mean_molecular_weight = source.mean_molecular_weight
        self.mix_diff_coeffs_mass = source.mix_diff_coeffs_mass

        # Rate properties
        self.heat_release_rate = source.heat_release_rate
        self.net_production_rates = source.net_production_rates
        self.heat_production_rates = source.heat_production_rates
        self.net_rates_of_progress = source.net_rates_of_progress

        # Computed properties
        self.Lt = self.thermal_flame_thickness()
        self.Ld = self.diffusive_flame_thickness()
        self.Z = self.local_Z()
        self.phi = self.local_phi()

    def thermal_flame_thickness(self):
        """Compute the thermal flame thickness"""
        Tu = self.T[0]
        Tb = self.T[-1]
        return (Tb - Tu) / np.max(np.gradient(self.T, self.grid))

    def diffusive_flame_thickness(self):
        """Compute the diffusive flame thickness"""
        kcond = self.thermal_conductivity[0]
        cp = self.cp_mass[0]
        rho = self.density[0]
        k_diff = kcond / (rho * cp)
        return k_diff / self.velocity[0]

    def global_progress(self):
        pass
        

    def local_Z(self):
        """Compute the local mixture fraction """
        Z = []
        for i, _ in enumerate(self.grid):
            self.gas.TPY = self.T[i], self.P, self.Y.T[i]
            Z_loc = self.gas.mixture_fraction(self.fuel, self.air)
            Z.append(Z_loc)
        return Z

    def local_phi(self):
        """Compute the local equivalence ratio"""
        phi = []
        for i, _ in enumerate(self.grid):
            self.gas.TPY = self.T[i], self.P, self.Y.T[i]
            phi_loc = self.gas.equivalence_ratio(self.fuel, self.air)
            phi.append(phi_loc)
        return phi
