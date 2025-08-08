import matplotlib.pyplot as plt
import numpy as np

""" A simple plot lib. Users may need to adjust final results."""


def plot_result(
    data,
    xlabel,
    ylabel,
    figtitle,
    nbins=30,
    limits=None,
    density=False,
    savepath=None,
    cumulative=None,
):
    if (isinstance(data, np.ndarray) and len(data.shape) > 1):
        print("[ERROR]: input data should be one dimensional. ")
        return

    f = plt.figure()
    ax = plt.gca()
    params = {
        "legend.fontsize": 22,
        "axes.labelsize": 22,
        "axes.titlesize": 22,
        "xtick.labelsize": 22,
        "ytick.labelsize": 22,
    }
    plt.ticklabel_format(style="sci", axis="x", scilimits=(0, 3))

    histogram = plt.hist(data, nbins, range=limits, density=density, cumulative=cumulative)

    ax.set_xlabel(xlabel, fontsize=22)
    ax.set_ylabel(ylabel, fontsize=22)
    ax.tick_params(axis="x", labelsize=18)
    ax.tick_params(axis="y", labelsize=18)
    ax.set_title(figtitle)

    if savepath is not None:
        plt.savefig(savepath, bbox_inches="tight")
    plt.show()
    
    return histogram
