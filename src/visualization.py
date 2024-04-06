import os
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import MultipleLocator, AutoMinorLocator

def plt_metric(metrics, output):
    """
    Parameters
    -----------
    metrics : list 
    output : Path object

    """
    #plt.rc('text', usetex=True)
    plt.rc('font', family='serif')
    plt.rc('xtick',labelsize=26)
    plt.rc('ytick',labelsize=26)
    plt.rc('grid', c='0.5', ls='-', alpha=0.5, lw=0.5)
    fig = plt.figure(figsize=(20,16))
    ax = fig.add_subplot(1,1,1)
    
    border_width = 1.5; axis_fs=44
    ax.spines['top'].set_linewidth(border_width)
    ax.spines['right'].set_linewidth(border_width)
    ax.spines['bottom'].set_linewidth(border_width)
    ax.spines['left'].set_linewidth(border_width)
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    ax.tick_params(direction='in', length=8, width=border_width, which='major', top=True, right=True)
    ax.tick_params(direction='in', length=4, width=border_width, which='minor', top=True, right=True)

    ax.set_xlabel(r'Graph Diameter Length, d', fontsize=axis_fs)
    ax.set_ylabel(r'Probability', fontsize=axis_fs)
    diams = []
    for metric in metrics:
        for m in metric: diams.append(m)

    weights = np.ones_like(diams)/float(len(diams))
    hist, edges, p = plt.hist(diams, weights=weights, bins=40, alpha=0.5, color='midnightblue',label='N1888+')
    #plt.legend()

    ax.set_xlim((0,max(diams)))
    plt.savefig(os.path.join(output,"diam.png"))
