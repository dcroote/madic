import os
from madic import command_line
import pytest
import subprocess


class TestCommandLine(object):
    # integration test

    def setup_method(self):
        self.parser = command_line.create_parser()

        self.dir_path = os.path.abspath(os.path.join(os.path.dirname(__file__),
                                                     '../../examples'))
        self.csv_path = os.path.join(self.dir_path,
                                     'madic_skyline_daily_data.csv')
        self.ref_path = os.path.join(self.dir_path,
                                     'madic_skyline_daily_reference.csv')

    def test_no_args(self):
        with pytest.raises(SystemExit):
            self.parser.parse_args()

    def test_command_line_execution(self, tmpdir):
        # from raw data to saved summary and data files

        out_summary_path = tmpdir.join('madic_summary.csv')
        out_log_path = tmpdir.join('madic.log')

        args = self.parser.parse_args([self.csv_path, self.ref_path, '_', '1',
                                       '--results_summary',
                                       '%s' % out_summary_path,
                                       '--log',
                                       '%s' % out_log_path])

        command_line.run_madic(args)

        with open(os.path.join(self.dir_path, 'madic_summary.csv')) as f:
            expected_summary = f.read()

        assert out_summary_path.read() == expected_summary

    def test_console_entry_point(self, tmpdir):

        out_summary_path = tmpdir.join('madic_summary.csv')
        out_log_path = tmpdir.join('madic.log')

        subprocess.call(['madic', self.csv_path, self.ref_path, '_', '1',
                         '--results_summary',
                         '%s' % out_summary_path,
                         '--log',
                         '%s' % out_log_path])

        with open(os.path.join(self.dir_path, 'madic_summary.csv')) as f:
            expected_summary = f.read()

        assert out_summary_path.read() == expected_summary
