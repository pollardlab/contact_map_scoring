
import itertools
from cooltools.lib.plotting import *


#### Helper Functions ####
def load_map(index, metadata, filepath):
    """
    Loads example and wt maps
    Input:
        index in metadata dataframe (ex. 0)
    Output:
        two numpy arrays (wt and variant)
    """
    example = metadata.iloc[index]
    wt_path = "{}{}_{}_{}_{}_wt.npy".format(filepath, example[0], example[1], example[2], example[3])
    del_path = "{}{}_{}_{}_{}_del.npy".format(filepath, example[0], example[1], example[2], example[3])
    wt_matrix = np.load(wt_path)
    del_matrix = np.load(del_path)
    title = "{}:{}-{}, {}".format(example[0], example[1],example[2],example[3])
    return wt_matrix, del_matrix, title



#### Plotting Functions ####
def pcolormesh_45deg(ax, mat, *args, **kwargs):
    #https://stackoverflow.com/questions/12848581/is-there-a-way-to-rotate-a-matplotlib-plot-by-45-degrees
    n = mat.shape[0]
    # create rotation/scaling matrix
    t = np.array([[1,0.5],[-1,0.5]])
    # create coordinate matrix and transform it
    A = np.dot(np.array([(i[1],i[0]) for i in itertools.product(range(n,-1,-1),range(0,n+1,1))]),t)
    # plot
    im = ax.pcolormesh(A[:,1].reshape(n+1,n+1),A[:,0].reshape(n+1,n+1),np.flipud(mat),*args, **kwargs)
    im.set_rasterized(True)
    ax.set_ylim(0,n)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.spines['left'].set_visible(False)
    _ = ax.set_xticks([])
    _ = ax.set_yticks([])
    ax.plot([0, n/2], [0, n], 'k-',linewidth=1)
    ax.plot([n/2, n], [n, 0], 'k-',linewidth=1)
    ax.set_aspect(.5)
    return im


def simple_plot(map_a, map_b, loops_a=None, loops_b=None, TAD_a=None, TAD_b=None, title = ""):
    """
    Simple plots of two maps
    Add tracks of insulation, points of loops
    """
    plt.rcParams['figure.facecolor'] = 'white'
    plt.rcParams['pdf.fonttype'] = 42
    plt.rcParams['font.size']= 12
    vmin, vmax= -2, 2
    plot_width, plot_width1D = 3, 0.75
    fig, gs = gridspec_inches([plot_width, .75, plot_width], [plot_width])
    #fig, gs = GridSpec([plot_width, .75, plot_width], [plot_width])

    # Plot Matrices
    plt.subplot(gs[0,0])
    im = plt.matshow(map_a, fignum=False, cmap= 'RdBu_r', vmax=vmax, vmin=vmin)
    plt.xticks([])
    if TAD_a is not None:
        for xc in TAD_a:
            plt.axvline(x=xc,color='green')
    if loops_a is not None:
        for c in loops_a:
            plt.scatter(c[1],c[0],color='red',s=10)

    plt.subplot(gs[0,2])
    im = plt.matshow(map_b, fignum=False, cmap= 'RdBu_r', vmax=vmax, vmin=vmin)
    plt.xticks([])
    if TAD_b is not None:
        for xc in TAD_b:
            plt.axvline(x=xc,color='green')
    if loops_b is not None:
        for c in loops_b:
            plt.scatter(c[1],c[0],color='red',s=10)

    plt.suptitle(title, y = 1.1)
    plt.show()