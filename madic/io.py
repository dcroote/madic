from __future__ import division
import pandas as pd
import numpy as np
import logging
from scipy.signal import savgol_filter
from madic import utils


logger = logging.getLogger(__name__)


def replicate_to_sample_name(series, delimiter, delimiter_pos):
    """ Extracts sample name from rep column based on delimiter and position

    Args:
        series (pd.Series): Replicate names
        delimiter (str): string by which to split replicate name e.g. '_'
        delimiter_pos (int): delimiter position that identifies a sample name

    Returns:
        pd.Series containing extracted sample names

    Example:
        >> rep = pd.Series(['Study1_Site4_ConditionA_injection_094'])
        >> print(replicate_to_sample_name(rep, '_', 2).values[0])
        'ConditionA'
    """

    result_series = series.str.split(delimiter, expand=True)[delimiter_pos]
    result_series.name = 'sample_name'

    if len(result_series.unique()) == 1:
        logger.warn('Warning: only one sample detected: did you specify '
                    '`delimiter` and `delimiter-pos` '
                    'when loading data?')

    return result_series


def expand_comma_sep_series(series, smooth=False,
                                  window_length=5, polyorder=1):
    """ Expands comma separated series into array

    Necessary to convert a time and chromatogram transition data within a
    Skyline report into a useable format.

    Args:
        series (pd.Series): Series of strings containing comma separated floats
        smooth (bool):  whether to perform smoothing (default=False)
        window_length (int): Savitzky-Golay filter window length (used only
            when smooth=True)
        polyorder (int): Savitzky-Golay filter polyorder (used only when
            smooth=True)

    Returns:
        pd.Series of np.arrays
    """

    expanded = series.apply(np.fromstring, sep=',')

    if smooth:
        expanded = expanded.apply(savgol_filter,
                                  window_length=window_length,
                                  polyorder=polyorder)

    return expanded


def read_transition_report(infile, delimiter=None, delimiter_pos=None,
                           smooth_chromatograms=True):
    """ Load Skyline transition report csv

    Note:
        csv file should be generated by exporting Skyline-daily data using the
        MADIC Skyline-daily report (located in:
        static/madic_skyline_daily_transition_report_and_chromatogram.skyr).
        This custom report can be loaded within Skyline as follows:
        File --> Export --> Report --> Edit list --> Import

        See :func:`io.replicate_to_sample_name` for `delimiter` and
        `delimiter_pos` use

    Args:
        infile (str): path to csv file
        delimiter (str): If not None, split replicate names by this/these
            character(s)
        delimiter_pos (int): refers to the list index of split replicate name
            array that identifies the sample name

    Returns:
        pd.DataFrame
    """

    logger.info("Reading transition report")

    df = pd.read_csv(infile)

    # more manageable column names
    df.rename(columns={'Replicate Name': 'rep',
                       'Acquired Time': 'acq_time',
                       'Protein Name': 'pro',
                       'Peptide Sequence': 'pep',
                       'Fragment Ion': 'prod_ion',
                       'Isotope Label Type': 'label',
                       'Retention Time': 'rt',
                       'Area': 'area',
                       'Start Time': 'rt_start',
                       'End Time': 'rt_end',
                       'Interpolated Times': 'times',
                       'Interpolated Intensities': 'intensities'
                       }, inplace=True)

    verify_transition_report_columns(df)

    # categorical dtype is more performant than object (~string) dtype
    df = utils.cols_to_category_dtype(df, ['pep', 'pro', 'prod_ion', 'label'])

    if delimiter is not None and delimiter_pos is not None:
        df['sample_name'] = replicate_to_sample_name(
            df.rep, delimiter, delimiter_pos)

    # add array columns derived from comma separated chromatogram columns
    df['intensities_arr'] = expand_comma_sep_series(df.intensities,
                                                    smooth=smooth_chromatograms)
    df['times_arr'] = expand_comma_sep_series(df.times)

    return df


def verify_transition_report_columns(df):
    """ Verify necessary columns are present

    Args:
        df: pd.DataFrame from transition report csv

    Raises:
        exception if any necessary columns are missing
    """
    necessary_columns = {'rep',
                         'acq_time',
                         'pro',
                         'pep',
                         'prod_ion',
                         'label',
                         'rt',
                         'area',
                         'rt_start',
                         'rt_end',
                         'times',
                         'intensities'}

    missing_cols = necessary_columns - set(df.columns)

    if len(missing_cols) != 0:
        raise ValueError(
            "Transition report missing columns: {}".format(missing_cols))


def write_out_data(df, outfile):
    """ Save data to csv

    Note:
        Repackages modified times and intensity arrays as comma separated
        strings

    Args:
        df: pd.DataFrame
    """

    logger.info("Writing all processed data to "
                "file: {}".format(outfile))

    df = df.copy()

    df['times_arr_str'] = df.times_arr.apply(
        lambda x: ','.join(['%s' % val for val in x]))

    df['intensities_arr_str'] = df.intensities_arr.apply(
        lambda x: ','.join(['%s' % val for val in x]))

    df.drop(['intensities_arr', 'times_arr'], axis=1, inplace=True)

    df.to_csv(outfile, index=False)
