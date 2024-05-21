#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# This code had been adapted from https://github.com/pinellolab/dictys
# Wang, L., Trasanidis, N., Wu, T. et al. Dictys: dynamic gene regulatory 
## network dissects developmental continuum with single-cell multiomics. 
### Nat Methods 20, 1368–1378 (2023).

import numpy as np

def auc(curves, times):
    """
	Computes area under the curves.
	Parameters
	----------
	times:	numpy.ndarray(shape=(n,))
		X coordinates
	curves:	numpy.ndarray(shape=(ny,n))
		Y coordinates, one for each y curve
	Returns
	-------
	numpy.ndarray(shape=(ny,))
		Area under the curves
	"""

    if len(times) < 2 or not (times[1:] > times[:-1]).all():
        raise ValueError('times must be increasing and have at least 2 values.')
    timesdiff = times[1:] - times[:-1]
    curvesmean = (curves[:,1:] + curves[:, :-1]) / 2
    ans = curvesmean@timesdiff
    return ans

def _curvesnamic_network_char_terminal_logfc_(curves, times):
	"""
	Computes terminal logFC for curves.
	Parameters
	----------
	times:	numpy.ndarray(shape=(n,))
		X coordinates
	curves:	numpy.ndarray(shape=(ny,n))
		Y coordinates, one for each y curve
	Returns
	-------
	numpy.ndarray(shape=(ny,))
		Terminal logFCs
	"""
	if len(times) < 2 or not (times[1:] > times[:-1]).all():
		raise ValueError('times must be increasing and have at least 2 values.')
	return (curves[:,-1] - curves[:,0]) / (times[-1] - times[0])

def _curvesnamic_network_char_transient_logfc_(curves, times):
    """
	Computes transient logFC for curves.
	Parameters
	----------
	times:	numpy.ndarray(shape=(n,))
		X coordinates 
	curves:	numpy.ndarray(shape=(ny,n))
		Y coordinates, one for each y curve
	Returns
	-------
	numpy.ndarray(shape=(ny,))
		Transient logFCs
    """
    n = curves.shape[1]
    times = (times - times[0]) / (times[-1] - times[0])
    curves = curves - np.median([curves, np.repeat(curves[:, [0]], n, axis = 1), np.repeat(curves[:, [-1]], n, axis = 1)], axis = 0)
    return auc(curves, times)

def _curvesnamic_network_char_switching_time_(curves, times):
    """
   	Computes switching time for curves.
   	Parameters
   	----------
   	times:	numpy.ndarray(shape=(n,))
   		X coordinates
   	curves:	numpy.ndarray(shape=(ny,n))
   		Y coordinates, one for each y curve
   	Returns
   	-------
   	numpy.ndarray(shape=(ny,))
   		Switching time
    """
    n = curves.shape[1]
    times = (times - times[0]) / (times[-1] - times[0])
    curves = np.median([curves, np.repeat(curves[:, [0]], n, axis = 1), np.repeat(curves[:, [-1]], n, axis = 1)], axis = 0)
    return (auc((curves.T - curves[:, -1]).T, times)) / (curves[:, 0] - curves[:, -1] + 1E-300)


def _curvesnamic_network_char_area_(curves, times):
    return abs(curves[:, 0] + curves[:, -1]) * abs(times[-1] - times[0])


