import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
import numpy as np


mpl.rcParams['axes.prop_cycle'] = mpl.cycler(
    color=["#C44DA5", '#4DC46D', "#4D69C4", "#E8AE5A", "#696969", "#F25E65"])


class _chromatogram_plot():

    def __init__(self, df, ax, rt_bounds, legend, rt_shift,
                 plot_kwargs, legend_kwargs, ignore_interference):

        # attributes
        self.df = df.copy()
        self.ax = ax
        self.rt_bounds = rt_bounds
        self.legend = legend
        self.rt_shift = rt_shift
        self.plot_kwargs = {} if plot_kwargs is None else plot_kwargs.copy()
        self.legend_kwargs = legend_kwargs
        self.ignore_interference = ignore_interference

    def add_legend(self):
        """ Adds a legend, with optional kwargs """

        # some legend defaults if no kwargs are passed
        if self.legend_kwargs is None:
            self.legend_kwargs = {'bbox_to_anchor': (1, 0.8),
                                  'loc': 1,
                                  'handlelength': 0.8}

        self.ax.legend(self.df.prod_ion.values, **self.legend_kwargs)

    def add_rt_bounds(self):
        """ Draw vertical lines for peak boundaries """

        start_val = self.df.rt_start.values[0]
        end_val = self.df.rt_end.values[0]

        if self.rt_shift:
            start_val -= self.rt_shift
            end_val -= self.rt_shift

        self.ax.axvline(start_val, ls='-', c='grey', lw=1.5, zorder=0)
        self.ax.axvline(end_val, ls='-', c='grey', lw=1.5, zorder=0)

    def adjust_axes(self):
        """ Adjusts ticks, labels, and spines """

        # reduce number of ticks
        self.ax.yaxis.set_major_locator(MaxNLocator(2))
        self.ax.xaxis.set_major_locator(MaxNLocator(2))

        # log format for tick labels
        self.ax.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))

        self.ax.set_xlabel('Time')
        self.ax.set_ylabel('Intensity')

        # despine for aesthetics
        self.ax.spines['right'].set_visible(False)
        self.ax.spines['top'].set_visible(False)

    def add_transitions(self, row):
        """ Plot transition """

        if 'interference' in row and row.interference and \
                self.ignore_interference:
            return

        t = np.array(row.times_arr)  # explicit copy
        if self.rt_shift:
            t -= self.rt_shift

        self.ax.plot(t, row.intensities_arr, **self.plot_kwargs)

    def sort_dataframe(self):
        """ Sorts dataframe in descending transition order """

        self.df['_prod_ion_len'] = self.df.prod_ion.str.len()
        self.df.sort_values(['_prod_ion_len', 'prod_ion'],
                            ascending=[False, False],
                            inplace=True)

    def plot(self):
        """ Draw plot """

        self.sort_dataframe()

        self.df.apply(self.add_transitions, axis=1)

        if self.legend:
            self.add_legend()

        if self.rt_bounds:
            self.add_rt_bounds()

        self.adjust_axes()

        return self.ax


def plot_chromatogram(df, ax=None, rt_bounds=True,
                      legend=False, rt_shift=None, plot_kwargs=None,
                      legend_kwargs=None, ignore_interference=False):
    """ Plot a chromatogram

    Args:
        df (pd.DataFrame): Dataframe containing time, intensity and other
            optional data
        ax (optional): Axes object to use for drawing the plot
        rt_bounds (bool, optional): If True, plot peak boundaries as vertical
            lines
        legend (bool, optional): If True, draw legend
        rt_shift (float, optional): Retention time offset, useful for
            correcting retention time shift if plotting chromatogram overlays
        plot_kwargs (dict, optional): Passed to ax.plot()
        legend_kwargs (dict, optional): Passed to ax.legend() if `legend=True`

    Returns:
        Axes object with drawing
    """

    if ax is None:
        f, ax = plt.subplots()

    chrom = _chromatogram_plot(df, ax, rt_bounds, legend, rt_shift,
                               plot_kwargs, legend_kwargs, ignore_interference)

    chrom.plot()

    return ax


def plot_panel_array(data, ref=None, **kwargs):
    """ Plot multiple subpanels of chromatogram data

    Peptides are plotted as rows and replicates as columns. If `ref` is
    provided, the first column of the figure will contain reference
    chromatograms

    Args:
        data (pd.DataFrame): chromatogram data
        ref (pd.DataFrame, optional): chromatogram data from standards or
            synthetic peptides to serve as a reference

    Returns:
        Matplotlib Figure object
    """
    data = data.copy()

    reps = sorted(data.rep.unique())
    peps = sorted(data.pep.unique())

    nrep = len(reps)
    npep = len(peps)

    if ref is not None:
        nrep += 1
        # arbitrary selection of first replicate from reference data
        selected_rep = ref.rep.unique()[0]

    f, ax0 = plt.subplots(npep, nrep, figsize=(2.5*nrep, 2.5*npep))
    ax = ax0.ravel()

    ind = 0
    for pep in peps:
        for rep in reps:

            # if reference is provided, plot in first column
            if ref is not None and ind % nrep == 0:
                subset = ref[(ref.pep == pep)
                             & (ref.rep == selected_rep)
                             & (ref.label == 'light')]
                plot_chromatogram(subset, ax=ax[ind], legend=True, **kwargs)
                ax[ind].set_title('Reference')
                ind += 1

            subset = data[(data.pep == pep)
                          & (data.rep == rep)
                          & (data.label == 'light')]
            plot_chromatogram(subset, ax=ax[ind], **kwargs)
            ax[ind].set_title(rep, y=1.1)

            if ind % nrep != 0:
                # remove y-axis labels in inner panels
                ax[ind].set_ylabel('')
            else:
                # add peptide name to y-axis label of first column
                ax[ind].set_ylabel(pep)
            ind += 1

    f.tight_layout()

    return f
