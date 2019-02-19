import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd
from madic import plotting
import numpy as np
from scipy import signal


class TestPlots(object):

    def setup_method(self):
        # dataframe with two peptides, 3 transitions each
        self.df = pd.DataFrame([['PEP1', 'y5'],
                                ['PEP1', 'y6'],
                                ['PEP1', 'y7'],
                                ['PEP2', 'y8'],
                                ['PEP2', 'y9'],
                                ['PEP2', 'y10']],
                               columns=['pep', 'prod_ion'])

        # co-eluting gaussians for both peptides (all 6 transitions)
        intensities = []
        for i in range(6):
            intensities.append(signal.gaussian(21, std=3)*(i+1))
        self.df['intensities_arr'] = intensities

        self.df['times_arr'] = [np.arange(21, dtype=np.dtype(float))]*6
        self.df['rt_start'] = 4
        self.df['rt_end'] = 16
        self.df['rep'] = 'rep1'
        self.df['label'] = 'light'

    def test_no_rt_shift(self):
        # confirm rt_shift does not change the underlying data

        single_pep = self.df[self.df.pep == 'PEP1'].copy()
        _ = plotting.plot_chromatogram(single_pep, rt_shift=0.1)

        original_time_arr = pd.Series([np.arange(21, dtype=np.float)]*6, name='times_arr')
        for ind, arr in single_pep.times_arr.iteritems():
            assert (arr == original_time_arr.loc[ind]).all()

    def test_plot_chromatogram(self, tmpdir):
        # plot chromatogram variants and save image

        single_pep = self.df[self.df.pep == 'PEP1'].copy()

        # test variants of the same chromatogram
        f, ax = plt.subplots()
        # no legend
        _ = plotting.plot_chromatogram(single_pep, ax=ax)
        # with RT shift
        _ = plotting.plot_chromatogram(single_pep, ax=ax, rt_shift=0.1)
        # with legend
        _ = plotting.plot_chromatogram(single_pep, ax=ax, legend=True)
        # with legend and kwargs
        legend_kwargs = {'bbox_to_anchor': (1, 0.5),
                         'loc': 1, 'title': 'my_legend'}
        _ = plotting.plot_chromatogram(single_pep, ax=ax, legend=True,
                                       legend_kwargs=legend_kwargs)

        path = tmpdir.join('test_single_chrom.pdf')

        f.savefig(str(path))

    def test_plot_panel_array(self, tmpdir):

        f = plotting.plot_panel_array(self.df)

        path = tmpdir.join('test_single_chrom.pdf')

        assert len(f.get_axes()) == 2
        f.savefig(str(path))
