import pandas as pd
from pandas.testing import assert_series_equal, assert_frame_equal
import numpy as np
from madic import qc, utils
from scipy import signal


class TestSignalToNoise(object):

    def setup_method(self):
        # 3 transitions as DataFrame rows
        # intensities vary by test
        self.df = pd.DataFrame({'rt_start': [3]*3,
                                'rt_end': [17]*3,
                                'times_arr': [np.arange(21)]*3})

    def test_calc_sn_simple_peak(self):
        # peak height (10) to median (1) = 10
        row = self.df.loc[0]
        intensities = [1, 1, 1, 1, 1, 1, 1, 1, 1, 5,
                       10,
                       5, 1, 1, 1, 1, 1, 1, 1, 1, 1]

        result = qc._calc_sn(row.times_arr,
                            np.array(intensities),
                            row.rt_start,
                            row.rt_end)

        assert result == 10

    def test_calc_sn_missing_rt_bounds(self):
        row = self.df.loc[0].copy()
        row['rt_start'] = np.nan
        row['rt_end'] = np.nan

        result = qc._calc_sn(row.times_arr,
                            np.arange(21),  # random values
                            row.rt_start,
                            row.rt_end)

        assert np.isnan(result)

    def test_calc_sn_zero_intensity(self):
        row = self.df.loc[0].copy()

        result = qc._calc_sn(row.times_arr,
                            np.zeros(21),
                            row.rt_start,
                            row.rt_end)

        assert np.isnan(result)

    def test_calc_sn_peak_but_zero_outside(self):
        row = self.df.loc[0].copy()
        intensities = np.zeros(21)
        intensities[10] = 5

        result = qc._calc_sn(row.times_arr,
                            intensities,
                            row.rt_start,
                            row.rt_end)

        assert np.isinf(result)

    def test_qc_sn_no_interference(self):
        # three co-eluting gaussian transitions with varying intensity
        intensities = []
        for i in range(3):
            intensities.append(signal.gaussian(21, std=3)*i)

        self.df['intensities_arr'] = intensities

        assert qc._sn_from_groupby(self.df, threshold_sn=3)

    def test_qc_sn_one_interference(self):
        # peptide is random noise, except
        # one transition has interference that should be ignored in the S/N
        intensities = []
        interference = []
        for i in range(3):
            if i == 1:
                intensities.append(signal.gaussian(21, std=3)*1000)
                interference.append(True)
            else:
                intensities.append(np.ones(21) + np.random.rand(21))
                interference.append(False)

        self.df['intensities_arr'] = intensities
        self.df['interference'] = interference

        assert ~qc._sn_from_groupby(self.df, threshold_sn=3)

    def test_qc_sn_all_interference(self):

        # all transitions have interference (edge case)
        intensities = []
        interference = []
        for i in range(3):
            intensities.append(signal.gaussian(21, std=3)*1000*i)
            interference.append(True)

        self.df['intensities_arr'] = intensities
        self.df['interference'] = interference

        assert ~qc._sn_from_groupby(self.df, threshold_sn=3)

    def test_qc_sn_missing_data(self):

        # no measured peptide intensity
        intensities = []
        for i in range(3):
            intensities.append(np.zeros(21))

        self.df['intensities_arr'] = intensities

        assert ~qc._sn_from_groupby(self.df, threshold_sn=3)

    def test_qc_sn_join(self):
        self.df['intensities_arr'] = [signal.gaussian(21, std=3)]*3
        self.df['rep'] = ['rep1' for _ in range(3)]
        self.df['pep'] = ['PEP1' for _ in range(3)]
        self.df['label'] = ['light', 'light', 'light']
        self.df['prod_ion'] = ['y8', 'y7', 'y6']

        categorical_cols = ['pep', 'prod_ion', 'label']
        for col in categorical_cols:
            self.df[col] = self.df[col].astype('category')

        expected_frame = self.df.copy()
        expected_frame['pass_signal_to_noise'] = True

        self.df = qc.eval_signal_to_noise(self.df)

        assert_frame_equal(self.df, expected_frame)


class TestTransitionRatio(object):

    def setup_method(self):
        self.df = pd.DataFrame([['rep1', 'PEP1', 'y5', 1000, 'light'],
                                ['rep1', 'PEP1', 'y6', 1000, 'light'],
                                ['rep1', 'PEP1', 'y7', 3000, 'light']],
                               columns=['rep', 'pep', 'prod_ion', 'area',
                                        'label'])

        self.expected_tr = [0.2, 0.2, 0.6]

        self.ref_df = self.df.copy()

        # ref ratio values added in tests
        self.ref_ratios = pd.DataFrame([['PEP1', 'y5', 'light'],
                                        ['PEP1', 'y6', 'light'],
                                        ['PEP1', 'y7', 'light']],
                                       columns=['pep', 'prod_ion', 'label'])

    def test_simple_column_normalization(self):
        # no groupby columns used

        result = utils.norm_col_sum_gb(self.df, 'area', grouping_cols=None)

        expected_series = pd.Series(self.expected_tr, name='area_frac')

        assert_series_equal(result, expected_series)

    def test_transition_ratio_heavy_and_light(self):
        # test transition ratio calculation without interference

        self.df['interference'] = False

        expected_series = pd.Series(self.expected_tr, name='area_frac')

        result = utils.norm_col_sum_gb(self.df,
                                       'area',
                                       ['rep', 'pep', 'label'])

        assert_series_equal(result, expected_series)

    def test_transition_ratio_with_interference(self):
        # test if single transition with interference is masked correctly

        self.df['interference'] = [False, False, True]

        expected_series = pd.Series([0.5, 0.5, np.nan], name='area_frac')

        result = utils.norm_col_sum_gb(self.df,
                                       'area',
                                       ['rep', 'pep', 'label'])

        assert_series_equal(result, expected_series)

    def test_transition_ratio_with_interference2(self):
        # test if single transition with interference is masked correctly

        self.df['interference'] = [False, False, True]

        expected_series = pd.Series([0.5, 0.5, np.nan], name='area_frac')

        result = utils.norm_col_sum_gb(self.df,
                                       'area',
                                       ['rep', 'pep', 'label'])
        assert_series_equal(result, expected_series)

    def test_transition_ratio_all_with_interference(self):
        # test masking if all transitions have interference

        self.df['interference'] = [True, True, True]

        expected_series = pd.Series([np.nan, np.nan, np.nan],
                                    name='area_frac')

        result = utils.norm_col_sum_gb(self.df,
                                            'area',
                                            ['rep', 'pep', 'label'])

        assert_series_equal(result, expected_series)

    def test_aggregate_reference_ratios(self):

        self.ref_ratios['ref_tr'] = self.expected_tr

        result = qc._aggregate_ref_data(self.ref_df)

        assert_frame_equal(result, self.ref_ratios)

    def test_qc_transition_ratio_absolute_within_limits(self):

        expected_frame = self.df.copy()
        expected_frame['transition_ratio'] = self.expected_tr
        expected_frame['ref_tr'] = self.expected_tr
        expected_frame['tr_individ_pass'] = [True, True, True]
        expected_frame['pass_transition_ratio'] = [True, True, True]

        dfout = qc.eval_transition_ratio(self.df, self.ref_df)

        assert_frame_equal(dfout, expected_frame)

    def test_qc_transition_ratio_absolute_outside_limits(self):
        self.df['area'] = [1000, 1000, 98000]

        expected_frame = self.df.copy()
        expected_frame['transition_ratio'] = [0.01, 0.01, 0.98]
        expected_frame['ref_tr'] = self.expected_tr
        expected_frame['tr_individ_pass'] = [False, False, False]
        expected_frame['pass_transition_ratio'] = [False, False, False]

        dfout = qc.eval_transition_ratio(self.df, self.ref_df)

        assert_frame_equal(dfout, expected_frame)

    def test_qc_transition_ratio_interference(self):
        # recompute ref transition ratio due to interference
        self.df['interference'] = [False, False, True]
        self.df['ref_tr'] = [0.2, 0.2, 0.6]

        result = utils.norm_col_sum_gb(self.df,
                                       'ref_tr',
                                       ['rep', 'pep', 'label'])

        expected_series = pd.Series([0.5, 0.5, np.nan],
                                    name='ref_tr_frac')

        assert_series_equal(result, expected_series)


class TestRecomputeRT(object):
    def setup_method(self):
        # two peaks with maxima at 6 and 5, respectively:
        intensities = [np.array([1.,  2., 3., 4., 5., 6., 5., 4., 3., 2., 1.]),
                       np.array([1., 2., 3., 4., 5., 4., 3., 2., 1.])]

        times = [np.arange(1, 12),
                 np.arange(1, 10)]

        self.df = pd.DataFrame({'intensities_arr': intensities,
                                'times_arr': times,
                                'rt_start': [4., 3.],
                                'rt_end': [8., 7.],
                                'rt': [4.1, 3.05],  # pretend wrong RT
                                })

    def test_recompute_correct_peak_rt(self):
        # testing to make sure retention times are corrected to 6 and 5

        expected_frame = self.df.copy()
        expected_frame['rt_original'] = expected_frame.rt.values
        expected_frame['rt'] = [6., 5.]

        result = qc.recompute_peak_rt(self.df)

        assert_frame_equal(result, expected_frame)

    def test_missing_peptide(self):
        self.df['rt_start'] = np.nan
        self.df['rt_end'] = np.nan
        self.df['rt'] = np.nan

        expected_frame = self.df.copy()
        expected_frame['rt_original'] = np.nan
        expected_frame['rt'] = np.nan

        result = qc.recompute_peak_rt(self.df)

        assert_frame_equal(result, expected_frame)

    def test_no_data_within_rt_bounds(self):
        self.df['rt_start'] = 6.001
        self.df['rt_end'] = 6.003
        self.df['rt'] = 6.002

        expected_frame = self.df.copy()
        expected_frame['rt_original'] = expected_frame.rt.values
        expected_frame['rt'] = np.nan

        result = qc.recompute_peak_rt(self.df)

        assert_frame_equal(result, expected_frame)


class TestRetentionTime():

    def setup_method(self):
        # heavy and light for two peptides
        self.df = pd.DataFrame([
            ['rep1', 'PEP1', 'y4', 5.51, 'light', False],
            ['rep1', 'PEP1', 'y5', 5.49, 'light', False],
            ['rep1', 'PEP1', 'y6', 5.51, 'light', False],
            ['rep1', 'PEP1', 'y4', 5.50, 'heavy', False],
            ['rep1', 'PEP1', 'y5', 5.50, 'heavy', False],
            ['rep1', 'PEP1', 'y6', 5.50, 'heavy', False],
            ['rep1', 'PEP2', 'y7', 9.51, 'light', False],
            ['rep1', 'PEP2', 'y8', 9.49, 'light', False],
            ['rep1', 'PEP2', 'y9', 9.51, 'light', False],
            ['rep1', 'PEP2', 'y7', 9.50, 'heavy', False],
            ['rep1', 'PEP2', 'y8', 9.50, 'heavy', False],
            ['rep1', 'PEP2', 'y9', 9.50, 'heavy', False]],
                               columns=['rep', 'pep', 'prod_ion', 'rt',
                                        'label', 'interference'])

    def test_pass_coelute(self):
        expected_frame = self.df.copy()
        expected_frame['pass_retention_time'] = True

        result = qc.eval_retention_time(self.df)

        assert_frame_equal(result, expected_frame)

    def test_transition_wrong_rt(self):
        # very wrong retention time for PEP1 y4, light peptide should fail
        self.df.loc[0, 'rt'] = 15

        expected_frame = self.df.copy()
        expected_frame['pass_retention_time'] = [False]*6 + [True]*6

        result = qc.eval_retention_time(self.df)

        assert_frame_equal(result, expected_frame)

    def test_transition_rt_with_interference(self):
        # despite very wrong retention time for PEP1 y4, should not fail
        # because it is identified as interference
        self.df.loc[0, 'rt'] = 15
        self.df.loc[0, 'interference'] = True

        expected_frame = self.df.copy()
        expected_frame['pass_retention_time'] = True

        result = qc.eval_retention_time(self.df)

        assert_frame_equal(result, expected_frame)

    def test_transition_rt_if_missing_one_heavy_peptide(self):
        # return np.nan for pass_retention_time for PEP2 because missing heavy
        self.df = self.df[~((self.df.pep == 'PEP2') &
                            (self.df.label == 'heavy'))]

        expected_frame = self.df.copy()
        expected_frame['pass_retention_time'] = [True]*6 + [False]*3

        result = qc.eval_retention_time(self.df)

        assert_frame_equal(result, expected_frame)

    def test_no_heavy_peptides(self):
        self.df = self.df[self.df.label == 'light']

        expected_frame = self.df.copy()
        expected_frame['pass_retention_time'] = False

        result = qc.eval_retention_time(self.df)

        assert_frame_equal(result, expected_frame)


class TestReplicate(object):
    def setup_method(self):
        # 2 replicate injections of a single peptide
        self.df = pd.DataFrame([
            ['sample1','rep1','PEP1','y5','light',True,True,True],
            ['sample1','rep1','PEP1','y6','light',True,True,True],
            ['sample1','rep1','PEP1','y7','light',True,True,True],
            ['sample1','rep2','PEP1','y5','light',True,True,True],
            ['sample1','rep2','PEP1','y6','light',True,True,True],
            ['sample1','rep2','PEP1','y7','light',True,True,True]
        ],
                               columns=['sample_name', 'rep', 'pep',
                                        'prod_ion', 'label',
                                        'pass_signal_to_noise',
                                        'pass_transition_ratio',
                                        'pass_retention_time'])

        categorical_cols = ['pep', 'prod_ion', 'label']
        for col in categorical_cols:
            self.df[col] = self.df[col].astype('category')

    def test_replicate_pass(self):
        expected_frame = self.df.copy()
        expected_frame['pass_all_replicate'] = True

        result = qc.eval_all_replicate(self.df)

        assert_frame_equal(result, expected_frame)

    def test_replicate_fail(self):
        # failing one filter should cause replicate QC to fail
        self.df.loc[0, 'pass_signal_to_noise'] = False

        expected_frame = self.df.copy()
        expected_frame['pass_all_replicate'] = False

        result = qc.eval_all_replicate(self.df)

        assert_frame_equal(result, expected_frame)


class TestSummarization(object):

    def setup_method(self):
        # two peptides, first passes all, second fails transition ratio
        # and retention time
        self.df = pd.DataFrame([
            ['sample1','rep1','PEP1','y5','light',True,True,True,True,False],
            ['sample1','rep1','PEP1','y6','light',True,True,True,True,False],
            ['sample1','rep1','PEP1','y7','light',True,True,True,True,False],
            ['sample1','rep2','PEP1','y5','light',True,True,True,True,False],
            ['sample1','rep2','PEP1','y6','light',True,True,True,True,False],
            ['sample1','rep2','PEP1','y7','light',True,True,True,True,False],

            ['sample1','rep1','PEP2','y5','light',True,True,False,True,False],
            ['sample1','rep1','PEP2','y6','light',True,True,True,True,False],
            ['sample1','rep1','PEP2','y7','light',True,True,True,True,False],
            ['sample1','rep2','PEP2','y5','light',True,False,True,True,False],
            ['sample1','rep2','PEP2','y6','light',True,False,True,True,False],
            ['sample1','rep2','PEP2','y7','light',True,False,True,True,False]
        ],
                               columns=['sample_name', 'rep', 'pep',
                                        'prod_ion', 'label',
                                        'pass_signal_to_noise',
                                        'pass_transition_ratio',
                                        'pass_retention_time',
                                        'pass_all_replicate',
                                        'interference'])

    def test_all_pass(self):
        expected_frame = pd.DataFrame([
            ['sample1','PEP1',True,True,True,True,False],
            ['sample1','PEP2',True,False,False,True,False]],
                                      columns=['sample_name', 'pep',
                                               'pass_signal_to_noise',
                                               'pass_transition_ratio',
                                               'pass_retention_time',
                                               'pass_all_replicate',
                                               'interference_corrected'])

        result = qc.summarize_results(self.df)

        assert_frame_equal(result, expected_frame)
