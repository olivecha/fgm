import pickle
from tqdm import tqdm
import multiprocessing
import numpy as np
import cantera as ct
from .pickle_flame import PickleFreeFlame

def solve_flame_at_phi(args):
    """
    Multiprocessing function solving one 1D flame
    """
    phi = args['phi']
    mech = args['mech']
    fuel = args['fuel']
    Tin = args['temp']
    P = args['P']
    air = args['air']
    gas = ct.Solution(mech)
    gas.TP = Tin, P
    gas.set_equivalence_ratio(phi, fuel, air)
    width = 0.16
    flame = ct.FreeFlame(gas, width=width)
    flame.set_refine_criteria(ratio=3.0, slope=0.1, curve=0.1)
    try:
        flame.solve(loglevel=0, auto=True)
        return PickleFreeFlame(flame, phi, fuel, air)
    except ct._utils.CanteraError:
        return

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
        pool = multiprocessing.Pool()
        inputs = []

        for phi in self.phi_values:
            fun_call = {'phi':phi,
                        'mech':self.mech,
                        'fuel':self.fuel,
                        'air':self.air,
                        'temp':self.T,
                        'P':self.P}
            inputs.append(fun_call)
        for f in tqdm(pool.imap(solve_flame_at_phi, inputs), total=ncalls):
            if f is not None:
                flames.append(f)
        self.flames = flames

    def save(self, filepath):
        """
        Save fgm to {filepath}
        """
        if len(filepath.split('.')) >= 2:
                pass
        else:
            filepath += '.fgm'
        with open(filepath, 'wb') as ffile:
            pickle.dump(self, ffile)




