import pandas as pd
from pandas.testing import assert_frame_equal
import numpy as np
from madic import interference
from scipy import signal


class TestInterference(object):

    def setup_method(self):
        # all pass_transition_ratio have been set to False for testing purposes
        self.df = pd.DataFrame({'rt_start': [3]*3,
                                'rt_end': [17]*3,
                                'times_arr': [np.arange(21)]*3,
                                'ref_tr': [0.2, 0.3, 0.4],
                                'pass_transition_ratio': [False, False, False]}
                               )

    def test_noise_should_not_be_interference(self):
        # random noise is not matrix interference
        intensities = [np.random.rand(21),
                       np.random.rand(21),
                       np.random.rand(21)]

        self.df['intensities_arr'] = intensities
        self.df['transition_ratio'] = [0.33, 0.33, 0.33]
        self.df['area'] = [3000, 3000, 3000]

        expected_frame = self.df.copy()
        expected_frame['interference'] = False

        result = interference.identify_interference(self.df)

        assert_frame_equal(result, expected_frame)

    def test_identify_interference(self):
        # 3 transitions, 1st has large interference
        intensities = [signal.gaussian(21, std=3)*1000,
                       signal.gaussian(21, std=3),
                       signal.gaussian(21, std=3)]

        self.df['intensities_arr'] = intensities
        self.df['transition_ratio'] = [0.98, 0.01, 0.01]
        self.df['area'] = [50000, 3000, 3000]

        expected_frame = self.df.copy()
        expected_frame['interference'] = [True, False, False]

        result = interference.identify_interference(self.df)

        assert_frame_equal(result, expected_frame)
