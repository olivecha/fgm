from tqdm import tqdm
import multiprocessing
import numpy as np
import cantera as ct
from .pickle_flame import PickleFreeFlame

def solve_flame_at_phi(args):
    """
    Multiprocessing function solving one 1D flame
    """
    phi = args[0]
    self = args[1]
    print(phi, self.fuel)
    gas = ct.Solution(self.mech)
    #gas.TPY = 300, 1*ct.one_atm, Y
    gas.TP = self.temp, self.P * ct.one_atm
    gas.set_equivalence_ratio(phi, self.fuel, self.air)
    width = 0.16
    flame = ct.FreeFlame(gas, width=width)
    flame.set_refine_criteria(ratio=3.0, slope=0.1, curve=0.1)
    flame.solve(loglevel=0, auto=True)
    return PickleFreeFlame(flame, self.fuel, self.air)


class FGMGenerator(object):
    """
    Class to generate FGM data
    by computing multiple FreeFlames in
    parallel
    """
    air = {"O2":0.21, "N2":0.79}
    width = 0.16

    def __init__(self, P=None, fuel=None, air=None, temp=None,
                 mech=None, phi_lo=None, phi_high=None, phi_step=None,
                 width=None):
        """
        Solves Cantera 1D flames to produce a FGM
        P: flame pressure (atm)
        fuel: dict with fuel molar fractions
        air: Specific air properties, defaults to {"O2":0.21, "N2":0.79}
        temp: Inlet temperature
        mech: path to kinetics file
        phi_lo: lower phi for parameter sweep
        """
        self.phi_values = np.arange(phi_lo, phi_high + phi_step/2, phi_step)
        self.P = P * ct.one_atm
        if air is not None:
            self.air = air
        if width is not None:
            self.width = width
        self.fuel = fuel
        self.T = temp
        self.gas = ct.Solution(mech)
        self.flames = None
        self.mech = mech
        self.temp=temp

        self.solve()

    def solve(self):
        """Solve Cantera FreeFlames in parallel"""
        flames = []
        ncalls = len(self.phi_values)
        with multiprocessing.Pool() as pool:
            for f in tqdm(pool.imap(solve_flame_at_phi, 
                                    zip(self.phi_values, [self]*ncalls)), total=ncalls):
                flames.append(f)
        self.flames = flames


        
        








