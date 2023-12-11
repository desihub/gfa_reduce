# The function fit_2dgaussian() used by GFA_reduce is removed from photutils v1.1.
# This script is basically a copy of the orignal fit_2dgaussian() from photutils v1.0.0 with
# minor changes to retain compatibility with newer versions of photutils (v1.1+), which is
# still needed by other scripts in the GFA_reduce package.

import warnings
from astropy.modeling import Fittable2DModel, Parameter
from astropy.modeling.fitting import LevMarLSQFitter
from astropy.modeling.models import Const2D, Gaussian2D
from astropy.utils.exceptions import AstropyUserWarning
import numpy as np
from photutils.morphology import data_properties

def fit_2dgaussian(data, error=None, mask=None):
    """
    Fit a 2D Gaussian plus a constant to a 2D image.

    Non-finite values (e.g., NaN or inf) in the ``data`` or ``error``
    arrays are automatically masked. These masks are combined.

    Parameters
    ----------
    data : array_like
        The 2D array of the image.

    error : array_like, optional
        The 2D array of the 1-sigma errors of the input ``data``.

    mask : array_like (bool), optional
        A boolean mask, with the same shape as ``data``, where a `True`
        value indicates the corresponding element of ``data`` is masked.

    Returns
    -------
    result : A `GaussianConst2D` model instance.
        The best-fitting Gaussian 2D model.
    """
    data = np.ma.asanyarray(data)

    if mask is not None and mask is not np.ma.nomask:
        mask = np.asanyarray(mask)
        if data.shape != mask.shape:
            raise ValueError('data and mask must have the same shape.')
        data.mask |= mask

    if np.any(~np.isfinite(data)):
        data = np.ma.masked_invalid(data)
        warnings.warn('Input data contains non-finite values (e.g., NaN or '
                      'infs) that were automatically masked.',
                      AstropyUserWarning)
    
    if error is not None:
        error = np.ma.masked_invalid(error)
        if data.shape != error.shape:
            raise ValueError('data and error must have the same shape.')
        data.mask |= error.mask
        weights = 1.0 / error.clip(min=1.e-30)
    else:
        weights = np.ones(data.shape)

    if np.ma.count(data) < 7:
        raise ValueError('Input data must have a least 7 unmasked values to '
                         'fit a 2D Gaussian plus a constant.')
    
    # assign zero weight to masked pixels
    if data.mask is not np.ma.nomask:
        weights[data.mask] = 0.

    mask = data.mask
    data.fill_value = 0.
    data = data.filled()
    
    # Subtract the minimum of the data as a rough background estimate.
    # This will also make the data values positive, preventing issues with
    # the moment estimation in data_properties. Moments from negative data
    # values can yield undefined Gaussian parameters, e.g., x/y_stddev.
    props = data_properties(data - np.min(data), mask=mask)

    init_const = 0.  # subtracted data minimum above
    init_amplitude = np.ptp(data)
    
    g_init = GaussianConst2D(constant=init_const, amplitude=init_amplitude,
                             x_mean=props.xcentroid,
                             y_mean=props.ycentroid,
                             x_stddev=props.semimajor_sigma.value,
                             y_stddev=props.semiminor_sigma.value,
                             theta=props.orientation.value)
    
    # original code from photutil.centroids.gaussian.py 1.0.0. Cannot be used in
    # its original form due to API changes in photutils v1.1
    #g_init = GaussianConst2D(constant=init_const, amplitude=init_amplitude,
    #                         x_mean=props.xcentroid.value,
    #                         y_mean=props.ycentroid.value,
    #                         x_stddev=props.semimajor_axis_sigma.value,
    #                         y_stddev=props.semiminor_axis_sigma.value,
    #                         theta=props.orientation.value)

    fitter = LevMarLSQFitter()
    y, x = np.indices(data.shape)
    gfit = fitter(g_init, x, y, data, weights=weights)

    return gfit


class GaussianConst2D(Fittable2DModel):
    """
    A model for a 2D Gaussian plus a constant.

    Parameters
    ----------
    constant : float
        Value of the constant.

    amplitude : float
        Amplitude of the Gaussian.

    x_mean : float
        Mean of the Gaussian in x.

    y_mean : float
        Mean of the Gaussian in y.

    x_stddev : float
        Standard deviation of the Gaussian in x. ``x_stddev`` and
        ``y_stddev`` must be specified unless a covariance matrix
        (``cov_matrix``) is input.

    y_stddev : float
        Standard deviation of the Gaussian in y. ``x_stddev`` and
        ``y_stddev`` must be specified unless a covariance matrix
        (``cov_matrix``) is input.

    theta : float, optional
        Rotation angle in radians. The rotation angle increases
        counterclockwise.
    """

    constant = Parameter(default=1)
    amplitude = Parameter(default=1)
    x_mean = Parameter(default=0)
    y_mean = Parameter(default=0)
    x_stddev = Parameter(default=1)
    y_stddev = Parameter(default=1)
    theta = Parameter(default=0)

    @staticmethod
    def evaluate(x, y, constant, amplitude, x_mean, y_mean, x_stddev,
                 y_stddev, theta):
        """Two dimensional Gaussian plus constant function."""

        model = Const2D(constant)(x, y) + Gaussian2D(amplitude, x_mean,
                                                     y_mean, x_stddev,
                                                     y_stddev, theta)(x, y)
        return model
