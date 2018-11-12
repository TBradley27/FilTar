import unittest
import sys
import pandas
from pandas.testing import assert_frame_equal

df = pandas.read_csv('/gpfs/afm/moxon/thomas2/APAtrap/validation/results/bam/sample/SRS1024051.bedgraph',
                       sep='\t',
                       names=['chromosome','start','stop','coverage'],
			dtype = {'chromosome': str, 'start': int, 'stop': int, 'coverage': float}
)

run1 = pandas.read_csv('/gpfs/afm/moxon/thomas2/APAtrap/validation/results/bam/run/SRR2146408.bedgraph',
                       sep='\t',
                       names=['chromosome','start','stop','coverage'],
                        dtype = {'chromosome': str, 'start': int, 'stop': int, 'coverage': float}
)

class testbedgraph(unittest.TestCase):

	def test_file_can_be_read(self):
		self.assertIsNotNone(df)
	def test_col_numbers_is_equal_to_4(self):
		col_nums = df.shape[1]
		self.assertEqual(col_nums,4)
	def test_row_numbers_is_equal_to_14438029(self):
		row_nums = df.shape[0]
		self.assertEqual(row_nums,14438029)
	def test_start_column_chrom1_is_monotonic(self):
		tmp_df = df.loc[df['chromosome'] == '1']
		bool = tmp_df['start'].is_monotonic
		self.assertTrue(bool)
	def test_stop_column_chrom1_is_monotonic(self):
		tmp_df = df.loc[df['chromosome'] == '1']
		bool = tmp_df['stop'].is_monotonic
		self.assertTrue(bool)
	def test_sample_coverage_equals_run1_coverage(self):
		#assert_frame_equal(df, run1, check_dtype=True)
		self.assertTrue(df.equals(run1))

if __name__ == '__main__':
	unittest.main()
