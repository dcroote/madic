from __future__ import division
import numpy as np
from madic import qc
import logging


logger = logging.getLogger(__name__)


def _req_tr_for_interference(vals):
    """ Scaling function for minimum transition ratio deemed interference

    Linear scaling between 0.7 and 0.95 for reference transition ratios
    between 0.4 and 0.7 with a floor of 0.7 and a ceiling of 0.95 otherwise

    This is useful because some peptides will naturally have skewed
    transition ratios due to a single strong transition. Therefore, for
    these transitions to qualify as interference, the threshold minimum
    transition ratio must be higher.

    Args:
        vals (int or np.array)

    Returns:
        Minimum transition ratio value(s) (int or np.array)
    """

    tr_ceil = 0.95
    tr_floor = 0.7
    # x_ceil = 0.7
    # x_floor = 0.4

    # precomputed y = mx + b
    # m = (tr_ceil - tr_floor) / (x_ceil - x_floor)
    # b = tr_floor - m*x_floor
    m = 0.833333
    b = 0.366667

    reg = vals*m + b

    return np.max([np.min([reg, tr_ceil]), tr_floor])


def identify_interference(df):
    """ Identify interference in MRM data

    Note:
        Operates under the following assumptions:
            * Interference will dominate transition ratio
            * Interference has greater intensity than general noise

    Args:
        df (pd.DataFrame): loaded using :func:`io.read_transition_report`

    Returns:
        Original pd.DataFrame with addition bool column indicating interference
    """

    df = df.copy()
    sub = df.copy()

    # limit to peptides failing transition ratio QC
    sub = sub[(~sub.pass_transition_ratio)]

    # avoid false positives of peptides with naturally unequal trans ratios
    # by scaling required trans ratio for interference according to the ref
    # trans ratio
    sub['min_tr_for_interference'] = sub.ref_tr.transform(
        _req_tr_for_interference)
    sub = sub[sub.transition_ratio > sub.min_tr_for_interference]

    # require distinguished peak (high signal to noise) for interference
    _tsn = []
    for _, row in sub.iterrows():
        _tsn.append(qc._calc_sn(row.times_arr,
                                row.intensities_arr,
                                row.rt_start,
                                row.rt_end))

    sub['transition_sn'] = _tsn
    sub = sub[sub.transition_sn > 10]

    # low threshold avoids false positive low intensity spikes
    # in an otherwise low intensity chromatogram
    sub = sub[sub.area > 2000]
    logger.info("Found {} instances of transition "
                "interference.".format(sub.shape[0]))

    # merge back into original dataframe
    df['interference'] = False
    df.loc[sub.index, 'interference'] = True

    return df
