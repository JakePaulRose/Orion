import itertools
import numpy as np
import pnumpy
from scipy.stats import median_absolute_deviation

def segment_times(timeseries, max_gap):
    """
    Returns an N-D array where each row represents a separate segmentation of continuous data with no gaps
    greater than the max gap.
    """
    time_segments = []
    is_contiguous = False
    arr_n = -1
    for i, t in enumerate(timeseries):
        if not is_contiguous:
            time_segments.append([t])
            arr_n += 1
        else:
            time_segments[arr_n].append(t)
        if i + 1 < len(timeseries):
            is_contiguous = (timeseries[i + 1] - t) < max_gap
    return time_segments

def match_dimensions(ndarray, onedarray):
    """
    Return an N-D array of shape ndarray.shape with the values of onedarray
    """
    ndreturn = []
    idx = 0
    for sublist in ndarray:
        num = len(sublist)
        subreturn = onedarray[idx : idx + num]
        idx += num
        ndreturn.append(subreturn)
    return ndreturn

def regularize(timeseries, dt, bin_about=None):
    """
    :param timeseries: defines max range of output
    :param dt: gap between outputted time points
    :param bin_about: define centre of bins. If outside of range of timeseries will extrapolate
    :return: regular timeseries with spacing dt defined within the range of input
    """
    if bin_about:
        return np.concatenate(
            [
                np.r_[bin_about : timeseries.min() : -dt][::-1],
                np.r_[bin_about + dt : timeseries.max() : dt],
            ]
        )
    else:
        return np.r_[timeseries.min() : timeseries.max() : dt]

def medsig(
    a: np.ndarray, include_zeros: bool = True, axis: int = None
):
    """
    Compute median and MAD-estimated scatter of array a
    :param a: numpy NDArray
    :param include_zeros: bool, default True. If False ignore zero values from median calculation
    :param axis: (int), default None. axis over which to calculate the median.
    """
    a = a if include_zeros else a[np.nonzero(a)]
    med = np.nanmedian(a, axis=axis)
    sig = median_absolute_deviation(a, axis=axis, nan_policy="omit")
    return med, sig

def rebin_err(t, f, dt=0.02, get_err_on_mean=False, bin_about=None):
    """
    Rebin a time-series with errors on the data (y-points).
    Apply unweighted average: ignore errors on the data to be binned and
                              perform a simple MAD estimation of the error
                              from the scatter of points
    """
    treg = regularize(t, dt=dt, bin_about=bin_about)
    nreg = len(treg)
    freg = np.zeros(nreg) + np.nan
    freg_err = np.zeros(nreg) + np.nan
    for i in np.arange(nreg):
        l = (t >= treg[i]) * (t < treg[i] + dt)
        if l.any():
            treg[i] = np.nanmean(t[l])
            freg[i], freg_err[i] = medsig(f[l])
            if get_err_on_mean:
                freg_err[i] /= np.sqrt(float(len(f[l])))
    l = np.isfinite(freg)
    return treg[l], freg[l], freg_err[l]

def rebin_err_chunks(t, f, dt, max_gap=0.02, get_err_on_mean=False, bin_about=0):
    """
    Re-bin a time series with errors on the data, but by chunking up into slices to avoid empty data.
    """
    times = segment_times(t, max_gap)
    fluxes = match_dimensions(times, f)
    times_binned = []
    fluxes_binned = []
    flux_errs_binned = []
    for i, each_time in enumerate(times):
        tbin, fbin, ferr = rebin_err(
            np.array(each_time),
            np.array(fluxes[i]),
            dt=dt,
            get_err_on_mean=get_err_on_mean,
            bin_about=bin_about,
        )
        times_binned.append(tbin)
        fluxes_binned.append(fbin)
        flux_errs_binned.append(ferr)
    treg = list(itertools.chain(*times_binned))
    freg = list(itertools.chain(*fluxes_binned))
    freg_err = list(itertools.chain(*flux_errs_binned))
    return treg, freg, freg_err