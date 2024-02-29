import ci_reduce.dark_current as dark_current
import matplotlib.pyplot as plt
import numpy as np

def plot_fit_and_measured():
    # measured values from DESI-3358
    t_celsius, e_per_pix_per_sec = dark_current.desi_3358_measurements()

    pad = 5
    t_celsius_grid = np.arange(np.min(t_celsius)-pad, np.max(t_celsius)+pad)

    dark_current_pred = dark_current.dark_current_rate(t_celsius_grid)

    plt.scatter(t_celsius, np.log(e_per_pix_per_sec), marker='+', s=100, 
                label='measured')
    plt.plot(t_celsius_grid, np.log(dark_current_pred), c='r', label='best fit')

    plt.xlabel('T (deg C)')
    plt.ylabel('ln(e-/pix/sec)')

    ax = plt.gca()
    ax.legend()
    plt.show()
