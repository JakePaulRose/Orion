import numpy as np
import matplotlib.pyplot as plt
from numpy.core.shape_base import _stack_dispatcher
from astropy.io import fits
import time
from astropy.stats import sigma_clip
from transitleastsquares import (
    transitleastsquares,
    cleaned_array,
    catalog_info,
    transit_mask
    )
import pandas as pd
from OrionPY import Tools

# Setting variables
min_mag = 9
max_mag = 15

# Load in data and obj_ids
hdu = Tools()
data = np.load('/data/jpr64/rebinned_data_test.npy', allow_pickle = True)
lc = [x for x in data]
obj_ids = [x for x in hdu.hdu[1].data['obj_id']]

# Calculating magnitudes
mean_mag = np.zeros(np.size(lc))
rms_mag = np.zeros(np.size(lc))

for i, j in enumerate(lc):
        magnitudes = -2.5 * np.log10(j[:,1]) + 20.2
        mean_mag[i] = np.nanmean(magnitudes)
        rms_mag[i] = np.sqrt(np.mean(np.square(magnitudes - mean_mag[i])))

# Keeping objects with rms_mag low enough    
index_keep = np.where((rms_mag > 10**(-2) & (mean_mag < max_mag) & (mean_mag > min_mag)))
index_keep = index_keep[0]
print(np.size(index_keep))

lc  = [lc[x] for x in index_keep]
obj_ids = [obj_ids[x] for x in index_keep]

# Looping over the remaining objects. and outputting them to a pandas df
output = pd.DataFrame()

for i, j in enumerate(lc):
    if i%50 == 0:
        output.to_pickle('/data/jpr64/TLS_search_undetrended.pkl')
    time, flux = cleaned_array(j[:,0], j[:,1]) # Clean of neg and nan values
    flux /= np.nanmean(flux) # Normalising
    sigma_clipped = sigma_clip(flux, sigma_lower = float('inf'), sigma_upper= 4) # Sigma clipping values to get rid of flares
    
    if np.size(time) < 100: continue # TLS fails if the array only has a few values in it. Anything less than 100 wouldn't be useful anyway.  
   
    # The actual TLS search
    model = transitleastsquares(time, flux)
    results = model.power(period_min = 1, period_max = 10, oversampling_factor = 2, use_threads = 40)
    
    # Adding the object id's to allow follow up
    results["obj_id"] = obj_ids[i]
    output = output.append(ignore_index=True) # Need to change this so it saves only the right columns

   
output.to_pickle('/data/jpr64/TLS_search_detrended.pkl')