import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import bin_tools
from OrionPY import OrionTools, Tools
import pandas as pd
from timeit import timeit

data = Tools()

flux_time_stack = np.stack((data.time, data.flux), axis = 2)

rebinned_data = []

# j = 0
start = timeit()
for i in flux_time_stack:
     rebinned_flux, rebinned_time, junk = bin_tools.rebin_err_chunks(i[:,0], i[:,1], dt=(1/240), max_gap=0.5)
     rebinned_data.append(np.stack((rebinned_time, rebinned_flux), axis = 1))
     
    #  j += 1
    #  if j == 10: break
end = timeit()

print(start-end)
x = np.asarray(rebinned_data, dtype = object)

"""Saving the file"""

np.save('/data/jpr64/rebinned_data.npy',x)

