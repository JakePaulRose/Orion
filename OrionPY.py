import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits

class OrionPy:
    def __init__(self, filename=None):
        if filename:
            self.hdu = self.load_file(filename)
        else:
            self.hdu = self.load_file('/data/jpr64/NG0535-0523_802_2017,2017S_CYCLE1807.fits')
     
    @staticmethod    
    def load_file(filename):
        return fits.open(filename)
   
    def load_flux(self, obj_index, calibrate_time = True):
        """Loads the time series and flux for a given obj_index (not id) then removes points with 0 flux. """
        obj_index = int(obj_index)
        self.time = self.hdu['hjd'].data[obj_index]
        self.time = self.time / (24*60*60)
        self.flux = self.hdu['sysrem_flux3'].data[obj_index]
        self.flux[self.flux == 0] = np.nan
    
        return self.time, self.flux
    
    def load_means(self):
        """Loads the mean flux, for every object in catalogue"""
        self.flux_mean = self.hdu['catalogue'].data['flux_mean']
        self.flux_rms = self.hdu['catalogue'].data['flux_wrms']
        self.mag_mean = self.hdu['catalogue'].data['mag_wmean']
        self.mag_rms = self.hdu['catalogue'].data['mag_rms']
        
        return self.flux_mean, self.mag_mean, self.mag_rmsgit 
    
    hello
    




