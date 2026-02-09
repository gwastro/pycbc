import os
import unittest
import tempfile
import numpy as np
from utils import simple_exit
from pycbc.io.hdf import HFile


class TestIOHDF(unittest.TestCase):

    def test_hfile_select_basic_and_premask(self):
        """Test HFile.select basic selection, premask as indices and boolean."""
        with tempfile.TemporaryDirectory() as td:
            p = os.path.join(td, "select.hdf")
            with HFile(p, 'w') as f:
                f.create_dataset('x', data=np.arange(10, dtype=np.int64))
                f.create_dataset('y', data=np.arange(10, dtype=np.int64) * 2)

            with HFile(p, 'r') as f:
                # simple select on x > 5
                idxs, (xs,) = f.select(lambda x: x > 5, 'x')
                np.testing.assert_array_equal(idxs, np.flatnonzero(np.arange(10) > 5))
                np.testing.assert_array_equal(xs, np.arange(6, 10))

                # premask as boolean array (only first 5 allowed)
                premask = np.zeros(10, dtype=bool)
                premask[:5] = True
                idxs2, _ = f.select(lambda x: x > 1, 'x', premask=premask)
                # only indices 2,3,4 should survive
                np.testing.assert_array_equal(idxs2, np.array([2, 3, 4], dtype=np.uint64))

                # premask as indices array
                premask_idx = np.array([7, 8, 9], dtype=int)
                idxs3, _ = f.select(lambda x: x > 7, 'x', premask=premask_idx)
                # only index 8,9 pass (x>7) while premask restricts to 7,8,9 -> final global indices 8 and 9
                np.testing.assert_array_equal(idxs3, np.array([8, 9], dtype=np.uint64))

    def test_hfile_select_mismatched_lengths_raises(self):
        """If datasets have different lengths, select should raise RuntimeError."""
        with tempfile.TemporaryDirectory() as td:
            p = os.path.join(td, "badlen.hdf")
            with HFile(p, 'w') as f:
                f.create_dataset('a', data=np.arange(5))
                f.create_dataset('b', data=np.arange(6))

            with HFile(p, 'r') as f:
                with self.assertRaises(RuntimeError):
                    f.select(lambda a, b: a > 0, 'a', 'b')

    def test_filedata_mask_and_get_column(self):
        """Test FileData.mask and get_column with a simple filter_func."""
        with tempfile.TemporaryDirectory() as td:
            p = os.path.join(td, "filedata.hdf")
            # create a file with a single top-level group so FileData can auto-select
            with HFile(p, 'w') as f:
                grp = f.create_group('grp')
                grp.create_dataset('a', data=np.arange(8))
                grp.create_dataset('b', data=np.arange(8) * 10)

            # Use the FileData class from the module under test
            from pycbc.io.hdf import FileData as FD
            fdata = FD(p)

            # Before setting filter_func, accessing mask should raise
            with self.assertRaises(RuntimeError):
                _ = fdata.mask

            # Now set a filter function that references 'a'
            fdata.filter_func = 'self.a > 4'
            # Access mask and column
            m = fdata.mask
            self.assertTrue(isinstance(m, np.ndarray) and m.dtype == bool)
            col = fdata.get_column('a')
            # Should return only values > 4
            np.testing.assert_array_equal(col, np.array([5, 6, 7]))

    def test_dictarray_save_and_reload(self):
        """Test that DictArray.save writes datasets and they can be reloaded."""
        from pycbc.io.hdf import DictArray
        with tempfile.TemporaryDirectory() as td:
            p = os.path.join(td, 'dictarray.hdf')
            data = {'a': np.array([1, 2, 3]), 'b': np.array([4, 5, 6])}
            da = DictArray(data=data)
            # ensure attrs exist to satisfy save implementation
            da.attrs = {'test': 'yes'}
            da.save(p)

            # open and verify datasets
            with HFile(p, 'r') as f:
                np.testing.assert_array_equal(f['a'][:], data['a'])
                np.testing.assert_array_equal(f['b'][:], data['b'])
                self.assertIn('test', f.attrs)

    def test_datafromfiles_get_column_concat(self):
        """Test DataFromFiles concatenates columns from multiple files."""
        from pycbc.io.hdf import DataFromFiles
        with tempfile.TemporaryDirectory() as td:
            p1 = os.path.join(td, 'f1.hdf')
            p2 = os.path.join(td, 'f2.hdf')

            # Create two files each with a single top-level group 'grp'
            with HFile(p1, 'w') as f:
                g = f.create_group('grp')
                g.create_dataset('val', data=np.array([1, 2, 3]))

            with HFile(p2, 'w') as f:
                g = f.create_group('grp')
                g.create_dataset('val', data=np.array([4, 5]))

            df = DataFromFiles([p1, p2], group='grp')
            out = df.get_column('val')
            np.testing.assert_array_equal(out, np.array([1, 2, 3, 4, 5]))

suite = unittest.TestSuite()
suite.addTest(unittest.TestLoader().loadTestsFromTestCase(TestIOHDF))


if __name__ == '__main__':
    results = unittest.TextTestRunner(verbosity=2).run(suite)
    simple_exit(results)
