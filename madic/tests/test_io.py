import os
import pandas as pd
import numpy as np
from pandas.testing import assert_series_equal
from madic import io


class TestChromatogramExpansion(object):

    def setup_method(self):
        # two rows of comma separated intensity chromatograms
        self.df = pd.DataFrame([['1,2,3,4,5,6,5,4,3,2,1'],
                                ['1,2,3,4,5,4,3,2,1']],
                               columns=['intensities'])

    def test_expand_comma_sep_series_no_smoothing(self):

        expected_series = pd.Series([np.array([1.,  2., 3., 4., 5., 6., 5.,
                                               4., 3., 2., 1.]),
                                     np.array([1., 2., 3., 4., 5., 4., 3.,
                                               2., 1.])],
                                    name='intensities')
        result = io._expand_comma_sep_series(self.df.intensities)

        assert_series_equal(result, expected_series)

    def test_expand_comma_sep_series_with_smoothing(self):

        expected_series = pd.Series([np.array([1.,  2., 3., 4., 4.6, 4.8,
                                               4.6, 4., 3., 2., 1.]),
                                     np.array([1., 2., 3., 3.6, 3.8, 3.6,
                                               3., 2., 1.])],
                                    name='intensities')

        result = io._expand_comma_sep_series(self.df.intensities,
                                                  smooth=True)

        assert_series_equal(result, expected_series)


class TestReplicateColumnSplit(object):

    def setup_method(self):
        self.series = pd.Series(['Site4_ConditionA_Part2_094',
                                 'Site4_ConditionA_Part3_095',
                                 'Site4_ConditionB_Part2_096',
                                 'Site4_ConditionB_Part3_097'
                                 ])

    def test_split_delimiter_position(self):

        expected_series = pd.Series(['ConditionA', 'ConditionA', 'ConditionB',
                                     'ConditionB'], name='sample_name')

        result = io.replicate_to_sample_name(self.series, '_', 1)

        assert_series_equal(result, expected_series)


def test_load_skyline_transition_report():

    report_path = os.path.abspath(os.path.join(os.path.dirname(__file__),
                                               '../../examples/'
                                               'madic_skyline_daily_data.csv'))
    df = io.read_transition_report(report_path,
                                   delimiter='_',
                                   delimiter_pos=1)

    assert sorted(df.label.unique()) == ['heavy', 'light']
    assert df.shape[0] == 40
    assert df.pep.unique().size == 2
    assert df.sample_name.unique().size == 2
    assert df.rep.unique().size == 4


def test_write_out_summary(tmpdir):
    summary = pd.DataFrame([
        ['sample1','PEP1',True,True,True,True,False],
        ['sample1','PEP2',True,False,False,True,False]],
                           columns=['sample_name', 'pep',
                                    'pass_signal_to_noise',
                                    'pass_transition_ratio',
                                    'pass_retention_time',
                                    'pass_all_replicate',
                                    'interference_corrected'])

    # write data
    path = tmpdir.join('summary.csv')
    summary.to_csv(str(path), index=False)

    # load expected file contents
    testsdir = os.path.abspath(os.path.dirname(__file__))
    expected_file = os.path.join(testsdir, 'data/for_testing_summary.csv')
    with open(expected_file) as f:
        expected = f.read()

    # compare contents
    assert path.read() == expected


def test_write_out_data(tmpdir):


    df = pd.DataFrame([
        ['rep1','PEP1','y5',True,True,True,True,False],
        ['rep1','PEP1','y6',True,True,True,True,False],
        ['rep1','PEP1','y7',True,True,True,True,False],
        ['rep2','PEP1','y5',True,True,True,True,False],
        ['rep2','PEP1','y6',True,True,True,True,False],
        ['rep2','PEP1','y7',True,True,True,True,False],

        ['rep1','PEP2','y5',True,True,False,True,False],
        ['rep1','PEP2','y6',True,True,True,True,False],
        ['rep1','PEP2','y7',True,True,True,True,False],
        ['rep2','PEP2','y5',True,False,True,True,False],
        ['rep2','PEP2','y6',True,False,True,True,False],
        ['rep2','PEP2','y7',True,False,True,True,False]
    ],
                      columns=['rep', 'pep', 'prod_ion',
                               'pass_signal_to_noise',
                               'pass_transition_ratio',
                               'pass_retention_time',
                               'pass_all_replicate',
                               'interference'])

    df['sample_name'] = 'sample1'
    df['label'] = 'light'
    df['times_arr'] = [np.arange(3)]*12
    df['intensities_arr'] = [[500.1, 800.9, 500.1]]*12

    path = tmpdir.join('data.csv')
    io.write_out_data(df, str(path))

    testsdir = os.path.abspath(os.path.dirname(__file__))
    expected_file = os.path.join(testsdir, 'data/for_testing_data.csv')
    with open(expected_file) as f:
        expected = f.read()

    assert path.read() == expected
