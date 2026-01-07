import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import binned_statistic
from matplotlib import cm
from matplotlib.colors import Normalize 
from scipy.interpolate import interpn
from scipy.optimize import curve_fit


def density_scatter(x , y, ax=None, sort=True, bins=20, axis_labels=['x','y'], **kwargs )   :
    """
    Scatter plot colored by point density from 2d histogram
    from: https://stackoverflow.com/questions/20105364/how-can-i-make-a-scatter-plot-colored-by-density-in-matplotlib
    Args:
        x, y - variables (np.array, pd.Series)
        ax - figure axis to plot
        sort - True plots densest points on top
        bins - number of equal-sized x-bins
        axis_labels - list of str
    Returns:
        ax
    """
    if ax is None :
        fig , ax = plt.subplots()
    data , x_e, y_e = np.histogram2d( x, y, bins=bins, density=True )
    z = interpn( ( 0.5*(x_e[1:] + x_e[:-1]) , 0.5*(y_e[1:]+y_e[:-1]) ) , data , np.vstack([x,y]).T , method = "splinef2d", bounds_error = False)

    #To be sure to plot all data
    z[np.where(np.isnan(z))] = 0.0

    # Sort the points by density, so that the densest points are plotted last
    if sort :
        idx = z.argsort()
        x, y, z = x[idx], y[idx], z[idx]

    ax.scatter( x, y, c=z, **kwargs )
    ax.set_xlabel(axis_labels[0])
    ax.set_ylabel(axis_labels[1])
    
    norm = Normalize(vmin = np.min(z), vmax = np.max(z))
    cbar = plt.colorbar(cm.ScalarMappable(norm = norm), ax=ax)
    cbar.ax.set_ylabel('point density')

    return ax


# read example data

data = pd.read_csv('testdata.csv', sep=';')
x = data['Vol'].values
y = data['BMfol'].values

#x = x[0:20000]
#y = y[0:20000]

# plot density scatter
fig, ax = plt.subplots(1,1)
density_scatter(x, y, ax, axis_labels=['Vol (m3ha-1)', 'foliage mass (ton ha-1)'], alpha=0.1)

# use scipy.stats.binned_statistics to provide median y for 6 x bins
yb, be, bn = binned_statistic(x, y, statistic='median', bins=60,range=[np.quantile(x, [0.01, 0.99])])

#  xbin centers
xb = 0.5*(be[0:-1] + be[1:])
# select only bins where data   
f = np.where((np.isfinite(yb)))
xb = xb[f]; yb = yb[f]

# plot binned data
ax.plot(xb, yb, 'ro', label='binned')

# fit non-linear model to data

def bef_nonlinear(x, a, b):
    #non-linear, fits well to binned mVMI data
    y =  a * (1 - np.exp(-b*x))
    return y

popt, pcov = curve_fit(bef_nonlinear, xb, yb) 

# evaluate fit at xx and add to figure
xx = np.linspace(min(xb), max(xb), 100)
yf = bef_nonlinear(xx, *popt)
ax.plot(xx, yf, 'k-', linewidth=2, label='fit')
ax.set_title('y=%.2E [1 - exp(-%.2E x)]' % (popt[0], popt[1]))

ax.legend()

plt.show()


###