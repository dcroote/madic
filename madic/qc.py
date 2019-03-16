from __future__ import division
import pandas as pd
import numpy as np
import logging
from madic import utils


logger = logging.getLogger(__name__)


def _aggregate_ref_data(df):
    """ Calculates transition ratios and takes the mean of any replicate
    injections of reference transition data

    Args:
        df (pd.DataFrame): loaded using :func:`madic.io.read_transition_report`

    Returns:
        Aggregated pd.DataFrame with reference transition ratios
        for each peptide
    """

    df = df.copy()

    gbcols = ['rep', 'pep', 'label']
    df['ref_tr'] = utils.norm_col_sum_gb(df, 'area', gbcols)

    agg = df.groupby(
        ['pep', 'prod_ion', 'label']).ref_tr.mean().reset_index()

    return agg


def eval_transition_ratio(df, ref_df, method='absolute', tolerance=0.15):
    """ Evaluate transition ratios against reference as a quality control step

    Args:
        df (pd.DataFrame): loaded using :func:`madic.io.read_transition_report`
            Contains sample data to evaluate.
        ref_df (pd.DataFrame): loaded using
            :func:`madic.io.read_transition_report`
            Contains synthetic standards / QC samples from which transition
        method (str): 'absolute' or 'relative'
        tolerance (float): fractional difference in transition ratios between
            reference value and measured

    Returns:
        Original pd.DataFrame with 3 additional columns:
            * (float) transition ratio
            * (float) reference transition ratio
            * (bool) pass fail transition ratio relative to reference
    """
    logger.info("Evaluating transition ratios")

    df = df.copy()

    df['transition_ratio'] = utils.norm_col_sum_gb(df, 'area',
                                                   ['rep', 'pep', 'label'])

    ref_processed = _aggregate_ref_data(ref_df)

    df = df.merge(ref_processed)

    if 'interference' in df.columns and df.interference.any():
        # recompute reference transition ratios masking
        # transitions with interference
        df['ref_tr'] = utils.norm_col_sum_gb(df, 'ref_tr',
                                             ['rep', 'pep', 'label'])
    if method == 'absolute':
        tr_lower_thresh = df.ref_tr - tolerance
        tr_upper_thresh = df.ref_tr + tolerance

    elif method == 'relative':
        tr_lower_thresh = df.ref_tr - tolerance * df.ref_tr
        tr_upper_thresh = df.ref_tr + tolerance * df.ref_tr

    else:
        raise ValueError(
            "Transition ratio method must be relative or absolute")

    # True if ratio between upper & lower bounds, or masked by np.nan
    df['tr_individ_pass'] = (((df.transition_ratio < tr_upper_thresh) &
                              (df.transition_ratio > tr_lower_thresh)) |
                             pd.isnull(df.transition_ratio))

    # all transitions for each peptide must pass
    df['pass_transition_ratio'] = df.groupby([
        'rep', 'pep', 'label']).tr_individ_pass.transform('all')

    return df


def _calc_sn(times_arr, intensities_arr, rt_start, rt_end):
    """ Signal to noise (background) calculation

    Calculate the ratio of peak height within peak boundaries to median
    intensity outside of peak boundaries after summing all arrays
    (chromatograms)

    Args:
        times_arr (np.array)
        intensities_arr (np.array)
        rt_start (float): peak boundary start time
        rt_end (float): peak boundary end time

    Returns:
        float: peak height to median background
    """

    if np.isnan(rt_start) or np.isnan(rt_end):
        # cannot calculate signal to noise without peak boundaries
        return np.nan

    peak_inds = (times_arr > rt_start) & (times_arr < rt_end)

    if not peak_inds.any():
        # edge case where peak boundaries chosen based on heavy
        # do not encompass light data
        return 0

    peak_height = intensities_arr[peak_inds].max()
    median_outside_peak = np.median(intensities_arr[~peak_inds])

    if median_outside_peak > 0:
        # expected case
        return peak_height / median_outside_peak
    elif median_outside_peak <= 0 and peak_height <= 0:
        # entire intensity array is zero
        return np.nan
    else:
        # peak but no background
        return np.inf


def _sn_from_groupby(g, threshold_sn):
    """ Compute the peak height to median background of the smoothed sum of
        transitions for an input DataFrame

    Note:
        Transitions with interference are ignored

    Args:
        df (pd.DataFrame): loaded using :func:`madic.io.read_transition_report`
        threshold_sn (float or int): Minimum peak height multiple of median
        background to pass QC

    Returns:
        bool: True if signal to noise is greater than `threshold_sn`
    """

    g = g.copy()

    rt_start = g.rt_start.values[0]
    rt_end = g.rt_end.values[0]

    if np.isnan([rt_start, rt_end]).any():
        logger.warn('Cannot compute signal to noise for peptide: {} in\
        replicate: {} due to missing peak boundaries'.format(g.pep.values[0],
                                                             g.rep.values[0]))
        return False

    # drop any interference
    if 'interference' in g.columns:
        g = g[~g.interference]

    if g.empty:
        # all transitions have interference
        return False

    # smooth each transition chromatogram and sum within group
    smoothed_sum = g.intensities_arr.values.sum()

    # times are the same for each transition within a peptide
    t_vals = g.times_arr.values
    t_mean = t_vals.mean()

    # find max within RT bounds
    peak_to_median = _calc_sn(t_mean, smoothed_sum, rt_start, rt_end)
    pass_fail = peak_to_median > threshold_sn

    return pass_fail


def eval_signal_to_noise(df, threshold_sn=3):
    """ Evaluate peptide signal to noise as a quality control step

    Args:
        df (pd.DataFrame): loaded using :func:`madic.io.read_transition_report`
        threshold_sn (int or float): minimum ratio of peak height to median
            intensity outside of peak boundaries in order to pass QC.
            Calculated using summed transition chromatograms

    Returns:
        Original pd.DataFrame with additional `pass_signal_to_noise` column
    """
    logger.info("Evaluating signal to noise")

    df = df.copy()

    if 'pass_signal_to_noise' in df.columns:
        df.drop('pass_signal_to_noise', axis=1, inplace=True)

    sn_gb_cols = ['rep', 'pep', 'label']

    gb = df.groupby(sn_gb_cols).apply(_sn_from_groupby,
                                      threshold_sn=threshold_sn)
    gb.name = 'pass_signal_to_noise'

    df = df.join(gb, on=sn_gb_cols)

    # can lose category dtype on merge
    df = utils.cols_to_category_dtype(df, ['pep', 'label'])

    return df


def recompute_peak_rt(df):
    """ Recalculate retention time peak based on chromatogram data

    Operates row-wise on a dataframe copy containing transition chromatograms.
    Useful for fixing the raw retention times chosen by Skyline after smoothing
    a noisy chromatogram

    Note:
        Renames existing retention time column to 'rt_original'

    Args:
        df (pd.DataFrame): loaded using :func:`madic.io.read_transition_report`

    Returns:
        Original pd.DataFrame with updated `rt` and new `rt_original` columns
    """
    logger.info("Recalculating retention times using smoothed chromatograms")

    df = df.copy()

    if 'rt' in df.columns and 'rt_original' not in df.columns:
        df['rt_original'] = df.rt

    # TODO: use apply instead of iterrows
    new_rts = np.zeros(df.shape[0])
    for c, (_, row) in enumerate(df.iterrows()):

        t = row.times_arr
        v = row.intensities_arr

        if np.isnan(row.rt_start) or np.isnan(row.rt_end):
            new_rts[c] = np.nan
            continue

        # find max within RT bounds
        peak_inds = (t > row.rt_start) & (t < row.rt_end)

        if not peak_inds.any():
            # edge case where peak boundaries do not encompass data
            new_rts[c] = np.nan
            continue

        peak_height_loc = t[peak_inds][np.where(
            v[peak_inds] == v[peak_inds].max())]

        if peak_height_loc.size == 1:
            new_rts[c] = peak_height_loc
        else:
            new_rts[c] = np.median(peak_height_loc)

    df['rt'] = new_rts

    return df


def eval_retention_time(df, rt_threshold=0.07):
    """ Evaluate transition retention time as a quality control step

    Light transition retention times are evaluated against mean peak retention
    time of heavy transitions within a given threshold

    Args:
        df (pd.DataFrame): loaded using :func:`madic.io.read_transition_report`
        rt_threshold (float): max allowable RT difference (min) between light
            transitions and heavy RT mean. Likely needs to be tuned based on
            chromatography
    Returns:
        Original pd.DataFrame with new bool column. True if light
        peptide transitions are within `rt_threshold` of heavy labeled peptide
        transitions
    """
    logger.info("Evaluating transition retention times")

    df = df.copy()

    if 'pass_retention_time' in df.columns:
        df.drop('pass_retention_time', axis=1, inplace=True)

    input_cols = df.columns.tolist()

    unique_light = df[df.label == 'light'].pep.unique()
    unique_heavy = df[df.label == 'heavy'].pep.unique()

    if unique_heavy.size == 0:
        logger.warn("No heavy peptides detected, cannot evaluate transition "
                    "retention times")
        df['pass_retention_time'] = False
        return df

    missing_heavy = np.setdiff1d(unique_light, unique_heavy)
    if missing_heavy.size > 0:
        logger.warn("Cannot evaluate transition retention times for the "
                    "following peptides due to missing heavy versions: "
                    "{}".format(missing_heavy))

    # mean RT of heavy transitions
    heavy_rt_means = df[df.label == 'heavy'].groupby([
        'rep', 'pep']).rt.mean().reset_index(name='rt_heavy')

    # broadcast heavy RT back on to heavy and light transitions
    df = df.merge(heavy_rt_means, how='left')  # 'left' in case heavy missing

    df['individ_rt_pass'] = (np.abs(df.rt_heavy - df.rt) <= rt_threshold)

    if 'interference' in df.columns:
        # ignore transitions with interference
        df['individ_rt_pass'] = df.individ_rt_pass | df.interference

    # all transitions for each peptide must pass
    df['pass_retention_time'] = df.groupby([
        'rep', 'pep']).individ_rt_pass.transform('all')

    return df[input_cols+['pass_retention_time']]


def eval_all_replicate(df, eval_metrics=None):
    """ Evaluate quality control filters across replicates

    Serves as the final quality control step where all `eval_metrics` must pass
    for all replicates

    Args:
        df (pd.DataFrame): loaded using :func:`madic.io.read_transition_report`
        eval_metrics (list, optional): previously evaluated metrics for which
            replicates must pass.
            Can include: `pass_signal_to_noise`, `pass_transition_ratio`, and
            `pass_retention_time`

    Returns:
        Original pd.DataFrame with additional bool `pass_all_replicate`
        column indicating peptide passes all `eval_metrics` for a sample
    """
    logger.info("Evaluating quality control metrics for all replicates")

    df = df.copy()

    if 'pass_all_replicate' in df.columns:
        df.drop('pass_all_replicate', axis=1, inplace=True)

    if eval_metrics is not None:
        eval_cols = []
        for qc in eval_metrics:
            eval_cols.append('pass_{}'.format(qc))
    else:
        # use default
        eval_cols = ['pass_signal_to_noise',
                     'pass_transition_ratio',
                     'pass_retention_time']

    # groupby columns
    eval_cols += ['sample_name', 'pep']

    agg = df[eval_cols].groupby(['sample_name', 'pep']).all().all(
        axis=1).reset_index(name='pass_all_replicate')

    df = df.merge(agg, how='left')

    # can lose category dtype on merge
    df = utils.cols_to_category_dtype(df, ['pep'])

    return df


def eval_data(df, ref_df, signal_to_noise=True, transition_ratio=True,
              peak_rt_recompute=True, retention_time=True, replicate=True):
    """ Convenience wrapper for calculating quality control metrics

    Args:
        df (pd.DataFrame): loaded using :func:`madic.io.read_transition_report`
            Contains sample data to evaluate.
        ref_df (pd.DataFrame): loaded using
            :func:`madic.io.read_transition_report`
            Contains synthetic standards / QC samples from which transition
            ratios are derived
        signal_to_noise (bool): If True, evaluate signal to noise
        transition_ratio (bool): If True, evaluate transition_ratios
        peak_rt_recompute (bool): If True, recompute the peak retention times
            of all transitions using smoothed chromatogram data (helpful with
            noisy data)
        retention_time (bool): If True, evaluate light transition retention
            against heavy peptide retention time within peak boundaries
        replicate (bool): If True, evaluate quality control metrics across
            replicates

    Returns:
        pd.DataFrame with additional columns corresponding to desired quality
        control filters and metrics
    """
    if signal_to_noise:
        df = eval_signal_to_noise(df)

    if transition_ratio:
        df = eval_transition_ratio(df, ref_df)

    if peak_rt_recompute:
        df = recompute_peak_rt(df)

    if retention_time:
        df = eval_retention_time(df)

    if replicate:
        df = eval_all_replicate(df)

    return df


def summarize_results(df, level='sample_name'):
    """ Summarize quality control results

    Args:
        df (pd.DataFrame): loaded using :func:`madic.io.read_transition_report`
        level (str): Choices: 'sample_name' or 'rep'
            Whether to return summary on a sample or replicate basis

    Returns:
        pd.DataFrame: DataFrame containing final pass / fail quality control
        results by peptide and sample
    """

    eval_cols = [x for x in df.columns if x.startswith('pass_')]

    summarized = df.groupby([
        level, 'pep'])[eval_cols].agg('all').reset_index()

    if 'interference' in df.columns:
        ic = df.groupby(['sample_name',
                         'pep']).interference.agg('any').reset_index()
        ic.rename(columns={'interference': 'interference_corrected'},
                  inplace=True)

        summarized = summarized.merge(ic)

    return summarized
