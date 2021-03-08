import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import bin_tools
from OrionPY import OrionTools, Tools

data = Tools()

flux_time_stack = np.stack((data.time, data.flux), axis = 2)

rebinned_data = []

j = 0

for i in flux_time_stack:
     rebinned_flux, rebinned_time, junk = bin_tools.rebin_err_chunks(i[:,0], i[:,1], dt=(1/240), max_gap=0.5)
     rebinned_data.append(list(zip(rebinned_time,rebinned_flux)))
     j += 1
     if j == 3: break
