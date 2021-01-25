import numpy as np
from scipy.interpolate import CubicSpline

def _atv_splinefwhm(_radii, _profile):

    sind = np.argsort(_radii)
    radii = _radii[sind]
    profile = _profile[sind]
    
    nrad = len(radii)

    _max = np.max(profile)

    if profile[0] != _max:
        print('Warning: Profile peak is off-center!')
        return np.nan

    fac = 400

    splrad = np.min(radii) + np.arange(nrad*fac + 1, dtype=float)*(np.max(radii)-np.min(radii))/(nrad*fac)

    nspl = len(splrad)

    spline = CubicSpline(radii, profile)

    splprof = spline(splrad)

    found = False
    max_splprof = np.max(splprof)

    _i = None
    for i in range(nspl):
        if splprof[i] < 0.5*max_splprof:
            found = True
            _i = i
            break

    if (not found) or (_i < 2):
        print('Warning: Unable to measure FWHM!')
        return np.nan

    fwhm = splrad[_i] + splrad[_i-1]

    return fwhm
    

    
