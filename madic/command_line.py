import argparse
import logging
from madic import qc, io, interference


def create_parser():
    """ Return parser """

    parser = argparse.ArgumentParser(description='MAtrix-Dependent \
                                     Interference Correction for \
                                     mass spectrometry data',
                                     usage="%(prog)s [options] csv ref_csv "
                                     "delimiter delimiter_pos")

    parser.add_argument('csv',
                        type=str,
                        help='Skyline data transition report csv')

    parser.add_argument('ref_csv',
                        type=str,
                        help='Skyline reference transition report csv')

    parser.add_argument('delimiter',
                        type=str,
                        help='Character(s) by which to separate replicate name'
                        ' to yield sample name. Used in conjunction with '
                        'sample-delimiter-positions')

    parser.add_argument('delimiter_pos',
                        type=int,
                        help='[int] position within the replicate name that '
                        'identifies a sample. Used in conjunction with '
                        'delimiter. Remember, python lists are zero-indexed!')

    parser.add_argument('--results_summary',
                        type=str,
                        default='madic_summary.csv',
                        help='Path to write summary results.'
                        '(Default: %(default)s)')

    parser.add_argument('--log',
                        type=str,
                        default='madic.log',
                        help='Path to write log file'
                        '(Default: %(default)s in the current directory)')

    parser.add_argument('--data_outfile',
                        type=str,
                        help='(Optional) path, which if specified, is where '
                        'all processed data will be saved')

    return parser


def run_madic(args):
    """ Execute madic pipeline

    Args:
        args: populated namespace from the ArgumentParser.parse_args() method

    """

    logging.basicConfig(filename=args.log,
                        filemode='w',
                        level=logging.DEBUG,
                        format="%(asctime)s %(name)-18.18s "
                        "[%(levelname)-8.8s]  %(message)s"
                        )

    logger = logging.getLogger(__name__)
    logger.info("Executing MADIC")
    logger.debug("Args %s" % args)
    print("Running MADIC")

    df = io.read_transition_report(args.csv,
                                   args.delimiter,
                                   args.delimiter_pos)

    ref_df = io.read_transition_report(args.ref_csv)

    processed = qc.eval_data(df, ref_df)

    processed = interference.identify_interference(processed)

    if processed.interference.any():
        logger.info("Reprocessing data to correct interference.")
        processed = qc.eval_data(processed, ref_df)

    # TODO: reorder columns? or have pre/post correction qc all pass?
    summary = qc.summarize_results(processed)
    logger.info("Writing summary results to "
                "file: {}".format(args.results_summary))
    summary.to_csv(args.results_summary, index=False)

    if args.data_outfile is not None:
        io.write_out_data(processed, args.data_outfile)

    logger.info("Completed successfully.")
    print("MADIC completed successfully. See madic.log for execution details.")


def main():
    """ Main function to parse and pass args for execution """

    parser = create_parser()
    args = parser.parse_args()

    run_madic(args)


if __name__ == '__main__':
    main()
