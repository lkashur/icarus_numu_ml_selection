import numpy as np
import uproot
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator

class PlotDescription:
    def __init__(self, cfg):
        plt.style.use('plot_style.mplstyle')
        for k, v in cfg.items():
            setattr(self, k, v)

def load_histograms(rf, key):
    """
    Loads the requested histograms from a ROOT file.

    Parameters
    ----------
    rf: ROOT file (uproot)
        The opened uproot ROOT file that contains the histogram.
    key: list[str]
        The names (key) of the histogram.

    Returns
    -------
    contents: list[np.array]
        The bin contents of the histogram normalized by the number of events.
    xedges: list[np.array]
        The x-bin edges of the histogram.
    yedges: list[np.array]
        The y-bin edges of the histogram (if 2D).
    centers: list[np.array]
        The bin centers of the histogram.
    """
    if isinstance(key, str) and len(rf[key].to_numpy()) > 2:
        contents, xedges, yedges = rf[key].to_numpy()
        return contents, xedges, yedges
    else:
        contents, xedges, centers = list(), list(), list()
        for k in key:
            c, e = rf[k].to_numpy()
            contents.append(c)
            xedges.append(e)
            centers.append((e[1:] + e[:-1]) / 2.0)
        return contents, xedges, centers

def plot_flow(rf, desc):
    """
    Plots a population flow chart according to the plot description
    using an input ROOT file.

    Parameters
    ----------
    rf: ROOT file (uproot)
        The opened uproot ROOT file that conains the histogram.
    desc: PlotDescription
        An instance of the PlotDescription class that describes the plot details.

    Returns
    -------
    None.
    """
    figure = plt.figure(figsize=(8,6))
    ax = figure.add_subplot()
    npops = len(desc.pops)
    bar_size = 1.0 / (npops+1)
    ylocs = np.arange(len(desc.cuts)) * (bar_size * (npops + 1))
    for pi, pop in enumerate(desc.pops):
        count = [load_histograms(rf, [f'sCount{desc.direction.upper()}_{pop}_{c}',])[0][0][0] for c in desc.cuts][::-1]
        ax.barh(ylocs + (pi + npops * -0.5)*bar_size, count, align='edge', height=bar_size, label=desc.labels[pi])
        for ci, c in enumerate(count):
            plt.text(1.2*10**desc.span[0], ci - (1.0-pi)*bar_size, f'{c:.2}', color='white', va="center", fontsize=18, weight='bold')

    ax.set(yticks=ylocs, yticklabels=desc.clabels[::-1], ylim=[0 - bar_size, len(desc.cuts)])
    ax.set_xscale('log')
    ax.set_xlim(10**desc.span[0], 10**desc.span[1])
    ax.set_ylim(-0.5,4.5)
    ax.set_xlabel('Interactions')
    ax.legend()
    figure.suptitle(desc.title)
    figure.savefig(desc.save)

def plot_confusion(rf, desc):
    """
    Plots a confusion matrix according to the plot description
    using an input ROOT file.

    Parameters
    ----------
    rf: ROOT file (uproot)
        The opened uproot ROOT file that conains the histogram.
    desc: PlotDescription
        An instance of the PlotDescription class that describes the plot details.

    Returns
    -------
    None.
    """
    figure = plt.figure(figsize=(8,6))
    ax = figure.add_subplot()

    contents, xedges, yedges = load_histograms(rf, desc.var)
    contents /= np.sum(contents, axis=0)    
    x, y = np.meshgrid(xedges, yedges)
    pc = ax.pcolormesh(x, y, contents, cmap='Blues', vmin=0, vmax=1.0)
    cb = figure.colorbar(pc, ax=ax, ticks=0.25*np.arange(5), label=desc.clabel)
    cb.ax.set_yticklabels([f'{i:.0%}' for i in cb.get_ticks()])

    ax.set_aspect('equal')
    ax.set_xticks(np.arange(contents.shape[0])+0.5, desc.entries)
    ax.set_yticks(np.arange(contents.shape[1])+0.5, desc.entries, rotation=90, va='center')
    for i in range(len(xedges)-1):
        for j in range(len(yedges)-1):
            c = 'white' if contents[i,j] > 0.5 else 'black'
            ax.text(j+0.5, i+0.5, f'{contents[i, j]:.2%}', ha='center', va='center', color=c, size=18)
    ax.set_xlabel(desc.xlabel)
    ax.set_ylabel(desc.ylabel)
    figure.suptitle(desc.title)
    figure.savefig(desc.save)

def plot_histogram_1d(rf, desc):
    """
    Plots a 1D histogram according to the plot description
    using an input ROOT file.

    Parameters
    ----------
    rf: uproot.reading.ReadOnlyDirectory
        The uproot handle for the ROOT file containing the input
        spectra/histograms.
    desc: PlotDescription
        An instance of the PlotDescription class that describes the
        plot details.

    Returns
    -------
    None.
    """
    figure = plt.figure(figsize=(8,6))
    ax = figure.add_subplot()

    contents, edges, centers = load_histograms(rf, desc.var)
    ax.hist(centers, weights=contents, range=(edges[0][0], edges[0][-1]), bins=edges[0], label=desc.var, **desc.plot_kwargs)
    ax.set_xlim(edges[0][0], edges[0][-1])
    ax.set_xlabel(desc.xlabel)
    ax.set_ylabel(desc.ylabel)
    ax.xaxis.set_major_locator(MaxNLocator(integer=True))
    ax.legend()
    figure.suptitle(desc.title)
    figure.savefig(desc.save)
