import matplotlib.pyplot as plt
import matplotlib as mpl


params = {'legend.fontsize': 20,
            'axes.labelsize': 22,
            'axes.titlesize':22,
            'xtick.labelsize':22,
            'ytick.labelsize':22}  # support latex in labels
mpl.rcParams.update(params)
plt.rc('text', usetex=True)
plt.rc('text.latex', preamble=r'\usepackage{amsmath}  \usepackage{amssymb} ')

""" A simple plot lib. Users may need to adjust final results."""

markers_def = ["o", "p", "^", "x", "d"]
linestyles_def = [(0, ()), (0, (1, 1)), (0, (5, 5)), (0, (3, 5, 1, 5)), (0, (5, 1))]
colors_def = ["b", "r", "c", "m", "orange", "k", "darkgreen", "darkorchid"]

def plot_result(data, line_names, xlabel, ylabel, figtitle, savepath=None,
    logx=False, logy=False, limx=None, limy=None, legend_pos="best", grid=False,
    colors=colors_def, markers=markers_def, linestyles=linestyles_def, linewidth=3, ax=None, showlegend=True,
    figsize=None,xticks=None,yticks=None,xtick_labels=None,ytick_labels=None):

    if line_names is None:
        line_names = [""] * len(data)

    f = plt.figure(figsize=figsize)
    if ax is None:
        ax = plt.gca()
    
    plt.ticklabel_format(style='sci', axis='x', scilimits=(0,3))

    if len(data) > 5:
        # we have too many lines, use auto line styles
        colors = [None] * len(data)
        if markers is None:
            markers = [None] * len(data)
        linestyles = [None] * len(data)

    for line_num, line_pair in enumerate(data):
        try:
            _ = line_pair[1]
        except:
            raise ValueError("data is wrongly formatted, please make sure it is a list of pairs."
                             "Each pair is the x and y coordinates.")
        ax.plot(line_pair[0], line_pair[1], color=colors[line_num], linestyle=linestyles[line_num],
                marker=markers[line_num], 
                markersize=8, linewidth=linewidth, label=line_names[line_num])
    
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
        
    if logx:
        ax.set_xscale("log")
    if logy:
        ax.set_yscale("log")
    if limx is not None:
        ax.set_xlim(limx)
    if limy is not None:
        ax.set_ylim(limy)
    if grid:
        ax.grid()
    

    if xticks is not None:
        print(xticks)
        ax.set_xticks(xticks)
        ax.set_xticklabels(xtick_labels)
    else:
        ax.tick_params(axis='x')
    if yticks is not None:
        ax.set_yticks(yticks)
        ax.set_yticklabels(ytick_labels)
    else:
        ax.tick_params(axis='y')

    ax.set_title(figtitle)
    if showlegend:
      ax.legend(loc=legend_pos)

    if savepath is not None:
        plt.savefig(savepath, bbox_inches="tight")
    # plt.show()
    return ax
