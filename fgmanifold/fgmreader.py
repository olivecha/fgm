import pickle
import numpy as np

class FGManifoldManager(object):

    Y_O2_IN = 0.22844261066374683
    Y_O2_OUT = 0.13601618257237072

    def __init__(self, filepath):
        """
        Read FGM data from {filepath} and initialize the FGM structure
        """
        obj = pickle.load(open(filepath, 'rb'))
        for ky in obj.__dict__:
            setattr(self, ky, obj.__dict__[ky])

        self.phi_values = [f.phi[0] for f in self.flames]


    def flame_at_phi(self, phi):
        """
        Finds the 1D flame closest to the phi value
        """
        return self.flames[np.argmin(np.abs(np.array(self.phi_values) - phi))]

    def Sl_at_phi(self, phi):
        """
        Return the nominal laminar flame speed at phi
        """
        return self.flame_at_phi(phi).velocity[0]

    def Lt_at_phi(self, phi):
        """
        Return the thermal flame thickness at phi
        """
        return self.flame_at_phi(phi).Lt


    def C(self, flame):
        Y_O2 = flame.Y[flame.gas.species_index('O2')]
        return (Y_O2 - self.Y_O2_IN) / (self.Y_O2_OUT - self.Y_O2_IN)
        
    