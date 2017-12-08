import matplotlib
matplotlib.use('Agg')
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

        self.df['times_arr'] = [np.arange(21)]*6
        self.df['rt_start'] = 4
        self.df['rt_end'] = 16
        self.df['rep'] = 'rep1'
        self.df['label'] = 'light'

    def test_plot_chromatogram(self, tmpdir):
        # plot a single chromatogram and save it

        single_pep = self.df[self.df.pep == 'PEP1'].copy()
        ax = plotting.plot_chromatogram(single_pep)

        path = tmpdir.join('test_single_chrom.pdf')

        f = ax.get_figure()
        f.savefig(str(path))

    def test_plot_panel_array(self, tmpdir):

        f = plotting.plot_panel_array(self.df)

        path = tmpdir.join('test_single_chrom.pdf')

        assert len(f.get_axes()) == 2
        f.savefig(str(path))
