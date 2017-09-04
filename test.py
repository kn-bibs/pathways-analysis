import unittest

import pandas as pd
import numpy as np

from utils import read_csv,fold_change


class Test_utils(unittest.TestCase):
    def test_read_csv(self):
        with self.assertRaises(TypeError):
            read_csv()
        with self.assertRaises(TypeError):
            read_csv('sample_expression.tsv')
        con,case = read_csv('sample_expression.tsv',[3,4,5])
        self.assertIsInstance(con,pd.DataFrame)
        self.assertIsInstance(case,pd.DataFrame)
        self.assertAlmostEqual(con['Normal_sample_1'].loc[['BAD']].values[0], 301.126,places=3)
        self.assertAlmostEqual(case['Tumour_sample_1'].loc[['BAD']].values[0], 368.033,places=3)

    def test_fold_change(self):
        with self.assertRaises(TypeError):
            fold_change()
        con, case = read_csv('sample_expression.tsv', [3, 4, 5])
        with self.assertRaises(TypeError):
            fold_change(case)
        f_cs = fold_change(case,con)
        self.assertIsInstance(f_cs,pd.DataFrame)
        t = f_cs['Tumour_sample_1'].loc[['BAD']].values[0]
        self.assertAlmostEqual(t,1.126,places=3)
        f_cs = fold_change(case, con,log2=True)
        t2 = f_cs['Tumour_sample_1'].loc[['BAD']].values[0]
        self.assertAlmostEqual(t2, np.log2(t))
        print(f_cs)


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
