import pickle
from tqdm import tqdm
import multiprocessing
import numpy as np
import cantera as ct
from .pickle_flame import PickleCounterflowTwinPremixedFlame
from .fgmreader import FGManifoldManager

def solve_flame_at_phi_uin(args):
    """
    Multiprocessing function solving one 1D flame
    """
    phi = args['phi']
    mech = args['mech']
    fuel = args['fuel']
    Tin = args['temp']
    P = args['P']
    air = args['air']
    Lx = args['Lx']
    u_in = args['u_in']

    gas = ct.Solution(mech)
    gas.TP = Tin, P
    gas.set_equivalence_ratio(phi, fuel, air)
    mass_flux = gas.density * u_in
    sflame = ct.CounterflowTwinPremixedFlame(gas, width=Lx)
    sflame.set_refine_criteria(ratio=2, slope=0.3, curve=0.3, prune=0.05)
    sflame.reactants.mdot = mass_flux
    try:
        sflame.solve(auto=True, loglevel=0)
        if np.max(sflame.T > 500):
            return PickleCounterflowTwinPremixedFlame(sflame, phi, fuel, air)
        else:
            return
    except ct._utils.CanteraError:
        return

class FGMGeneratorStrained(object):
    """
    Class to generate FGM data
    by computing multiple Strained falmes in
    parallel
    """
    air = {"O2":0.21, "N2":0.79}

    def __init__(self, P=None, fuel=None, air=None, temp=None,
                 mech=None, fgm_path=None, strain_data=None,
                 flames_per=None):
        """
        Solves Cantera 1D flames to produce a FGM
        P: flame pressure (atm)
        fuel: dict with fuel molar fractions
        air: Specific air properties, defaults to {"O2":0.21, "N2":0.79}
        temp: Inlet temperature
        mech: path to kinetics file
        phi_lo: lower phi for parameter sweep
        """
        # Load free flame data
        self.sdata = np.loadtxt(strain_data)
        self.fgm_ref = FGManifoldManager(fgm_path)

        self.P = P * ct.one_atm
        if air is not None:
            self.air = air
        self.fuel = fuel
        self.T = temp
        self.gas = ct.Solution(mech)
        self.flames = None
        self.mech = mech
        self.temp=temp
        self.solve(self.create_solve_calls())

    def create_solve_calls(self, n=None):
        """Compute a list of call for the MP fun"""
        mp_calls = []

        if n is None:
            n = self.sdata.shape[0]
            # For each value of phi
            for i, phi in enumerate(self.sdata[:, 0]):
                u_in_max = self.sdata[i, 1]
                Sl = self.fgm_ref.Sl_at_phi(phi)
                Lt = self.fgm_ref.Lt_at_phi(phi)
                for u_in in np.linspace(Sl, u_in_max, n):
                    call = {'phi':phi,
                            'mech':self.mech,
                            'fuel':self.fuel,
                            'temp':self.temp,
                            'P':self.P,
                            'air':self.air,
                            'Lx':5*Lt,
                            'u_in':u_in}
                    mp_calls.append(call)
        return mp_calls

    def solve(self, solve_calls):
        """Solve Cantera FreeFlames in parallel"""
        flames = []
        ncalls = len(solve_calls)
        pool = multiprocessing.Pool()

        for f in tqdm(pool.imap(solve_flame_at_phi_uin, solve_calls),
                      total=ncalls):
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


