import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import bin_tools

class OrionTools:
    def __init__(self, obj_index=None, filename=None):
        if filename:
            self.hdu = self.load_file(filename)
        else:
            self.hdu = self.load_file('/data/jpr64/NG0535-0523_802_2017,2017S_CYCLE1807.fits')
            if obj_index is not None:
                self.load_flux(obj_index)
            else:
                self.load_means()
     
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
        self.mag = -2.5*np.log10(self.flux)
    
    def load_means(self):
        """Loads the mean flux, for every object in catalogue"""

        self.flux_mean = self.hdu['catalogue'].data['flux_mean']
        self.flux_rms = self.hdu['catalogue'].data['flux_wrms']
        self.mag_mean = self.hdu['catalogue'].data['NGTS_gmag']
        self.mag_rms = self.hdu['catalogue'].data['mag_rms']
        self.gaia_gmag = self.hdu['catalogue'].data['gaia_gmag']
        self.gaia_gmag_err = self.hdu['catalogue'].data['gaia_gmag_err']

    def rebin(self):
        """Returns rebinned time and flux set to 6 min bins."""

        self.time_rebin, self.mag_rebin, junk = bin_tools.rebin_err_chunks(self.time, self.mag, dt=(1/240), max_gap=0.5)
        

    def rebinned_vals(self):
        """Calculates the new rebinned mean and rms"""

        self.calc_mag_mean = np.nanmean(self.mag_rebin)
        self.calc_mag_rms = np.sqrt((np.nanmean(np.square(self.mag_rebin))/np.size(self.mag_rebin)))

    
class Tools:
    
    def __init__(self):
        self.hdu = self.load_file('/data/jpr64/NG0535-0523_802_2017,2017S_CYCLE1807.fits')
        self.load_data_series()
        self.load_means()

    @staticmethod
    def load_file(filename):
            return fits.open(filename)
    
    def load_data_series(self, calibrate_time = True):
        """Loads the time series and flux for a given obj_index (not id) then removes points with 0 flux. """
        self.time = self.hdu['hjd'].data
        self.time = self.time / (24*60*60)
        self.flux = self.hdu['sysrem_flux3'].data
        self.flux[self.flux == 0] = np.nan
        self.mag = -2.5*np.log10(self.flux)

    def load_means(self):
        """Loads the mean flux, for every object in catalogue"""
        self.flux_mean = self.hdu['catalogue'].data['flux_mean']
        self.flux_rms = self.hdu['catalogue'].data['flux_wrms']
        self.imag_mean = self.hdu['catalogue'].data['NGTS_imag']
        self.gmag_mean = self.hdu['catalogue'].data['NGTS_gmag']
        self.mag_rms = self.hdu['catalogue'].data['mag_rms']
        self.gaia_gmag = self.hdu['catalogue'].data['gaia_gmag']
        self.gaia_gmag_err = self.hdu['catalogue'].data['gaia_gmag_err']
    
    def rebin(self):
        """Returns rebinned time and flux set to 6 min bins."""

        self.time_rebin, self.mag_rebin, junk = bin_tools.rebin_err_chunks(self.time, self.mag, dt=(1/240), max_gap=0.5)
        
    def rebinned_vals(self):
        """Calculates the new rebinned mean and rms"""

        self.calc_mag_mean = np.nanmean(self.mag_rebin)
        self.calc_mag_rms = np.sqrt((np.nanmean(np.square(self.mag_rebin))/np.size(self.mag_rebin)))     
