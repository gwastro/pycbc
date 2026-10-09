import os
import unittest
import tempfile
import numpy as np
from utils import simple_exit, parse_args_cpu_only
from pycbc.io.hdf import HFile

parse_args_cpu_only("io.hdf")


class TestIOHDF(unittest.TestCase):

    def test_hfile_select_basic_and_premask(self):
        """Test HFile.select basic selection, premask as indices/boolean."""
        with tempfile.TemporaryDirectory() as td:
            p = os.path.join(td, "select.hdf")
            with HFile(p, "w") as f:
                f.create_dataset("x", data=np.arange(10, dtype=np.int64))
                f.create_dataset("y", data=np.arange(10, dtype=np.int64) * 2)

            with HFile(p, "r") as f:
                # simple select on x > 5
                idxs, (xs,) = f.select(lambda x: x > 5, "x")
                np.testing.assert_array_equal(
                    idxs, np.flatnonzero(np.arange(10) > 5)
                )
                np.testing.assert_array_equal(xs, np.arange(6, 10))

                # premask as boolean array (only first 5 allowed)
                premask = np.zeros(10, dtype=bool)
                premask[:5] = True
                idxs2, _ = f.select(lambda x: x > 1, "x", premask=premask)
                # only indices 2,3,4 should survive
                np.testing.assert_array_equal(
                    idxs2, np.array([2, 3, 4], dtype=np.uint64)
                )

                # premask as indices array
                premask_idx = np.array([7, 8, 9], dtype=int)
                idxs3, _ = f.select(lambda x: x > 7, "x", premask=premask_idx)
                # only index 8,9 pass (x>7) while premask restricts to 7,8,9
                # -> final global indices 8 and 9
                np.testing.assert_array_equal(
                    idxs3, np.array([8, 9], dtype=np.uint64)
                )

    def test_hfile_select_mismatched_lengths_raises(self):
        """If datasets have different lengths, select should raise error."""
        with tempfile.TemporaryDirectory() as td:
            p = os.path.join(td, "badlen.hdf")
            with HFile(p, "w") as f:
                f.create_dataset("a", data=np.arange(5))
                f.create_dataset("b", data=np.arange(6))

            with HFile(p, "r") as f:
                with self.assertRaises(RuntimeError):
                    f.select(lambda a, b: a > 0, "a", "b")

    def test_filedata_mask_and_get_column(self):
        """Test FileData.mask and get_column with a simple filter_func."""
        with tempfile.TemporaryDirectory() as td:
            p = os.path.join(td, "filedata.hdf")
            # create file with single top-level group for FileData auto-select
            with HFile(p, "w") as f:
                grp = f.create_group("grp")
                grp.create_dataset("a", data=np.arange(8))
                grp.create_dataset("b", data=np.arange(8) * 10)

            # Use the FileData class from the module under test
            from pycbc.io.hdf import FileData as FD

            fdata = FD(p)

            # Before setting filter_func, accessing mask should raise
            with self.assertRaises(RuntimeError):
                _ = fdata.mask

            # Now set a filter function that references 'a'
            fdata.filter_func = "self.a > 4"
            # Access mask and column
            m = fdata.mask
            self.assertTrue(isinstance(m, np.ndarray) and m.dtype == bool)
            col = fdata.get_column("a")
            # Should return only values > 4
            np.testing.assert_array_equal(col, np.array([5, 6, 7]))

    def test_dictarray_save_and_reload(self):
        """Test DictArray.save writes datasets and they can be reloaded."""
        from pycbc.io.hdf import DictArray

        with tempfile.TemporaryDirectory() as td:
            p = os.path.join(td, "dictarray.hdf")
            data = {"a": np.array([1, 2, 3]), "b": np.array([4, 5, 6])}
            da = DictArray(data=data)
            # ensure attrs exist to satisfy save implementation
            da.attrs = {"test": "yes"}
            da.save(p)

            # open and verify datasets
            with HFile(p, "r") as f:
                np.testing.assert_array_equal(f["a"][:], data["a"])
                np.testing.assert_array_equal(f["b"][:], data["b"])
                self.assertIn("test", f.attrs)

    def test_datafromfiles_get_column_concat(self):
        """Test DataFromFiles concatenates columns from multiple files."""
        from pycbc.io.hdf import DataFromFiles

        with tempfile.TemporaryDirectory() as td:
            p1 = os.path.join(td, "f1.hdf")
            p2 = os.path.join(td, "f2.hdf")

            # Create two files each with a single top-level group 'grp'
            with HFile(p1, "w") as f:
                g = f.create_group("grp")
                g.create_dataset("val", data=np.array([1, 2, 3]))

            with HFile(p2, "w") as f:
                g = f.create_group("grp")
                g.create_dataset("val", data=np.array([4, 5]))

            dff = DataFromFiles([p1, p2], group="grp")
            res = dff.get_column("val")
            np.testing.assert_array_equal(res, np.array([1, 2, 3, 4, 5]))

    def test_dictarray_select_copy_vs_inplace(self):
        """Test DictArray.select with inplace=False and inplace=True."""
        from pycbc.io.hdf import DictArray

        # 1. Default copy (inplace=False)
        data = {
            "a": np.array([10, 20, 30, 40]),
            "b": np.array([100, 200, 300, 400])
        }
        da = DictArray(data=data)
        da_sub = da.select([1, 3])
        self.assertIsNot(da_sub, da)
        np.testing.assert_array_equal(da_sub.a, np.array([20, 40]))
        np.testing.assert_array_equal(da_sub.data["b"], np.array([200, 400]))
        # Original remains untouched
        self.assertEqual(len(da), 4)
        np.testing.assert_array_equal(da.a, np.array([10, 20, 30, 40]))

        # 2. Inplace mutation (inplace=True)
        ret = da.select([0, 2], inplace=True)
        self.assertIs(ret, da)
        self.assertEqual(len(da), 2)
        np.testing.assert_array_equal(da.a, np.array([10, 30]))
        np.testing.assert_array_equal(da.data["b"], np.array([100, 300]))

    def test_dictarray_remove_copy_vs_inplace(self):
        """Test DictArray.remove with copy, inplace, and empty index."""
        from pycbc.io.hdf import DictArray

        # 1. Default copy (inplace=False)
        da = DictArray(data={"a": np.array([1, 2, 3, 4, 5])})
        da_rem = da.remove([1, 3])
        self.assertIsNot(da_rem, da)
        np.testing.assert_array_equal(da_rem.a, np.array([1, 3, 5]))
        self.assertEqual(len(da), 5)

        # 2. Inplace mutation (inplace=True)
        ret = da.remove([0, 4], inplace=True)
        self.assertIs(ret, da)
        self.assertEqual(len(da), 3)
        np.testing.assert_array_equal(da.a, np.array([2, 3, 4]))

        # 3. Inplace remove with empty index
        ret_empty = da.remove([], inplace=True)
        self.assertIs(ret_empty, da)
        self.assertEqual(len(da), 3)

    def test_singledettriggers_apply_mask_scalar_and_empty(self):
        """Test SingleDetTriggers.apply_mask with scalar and empty indices."""
        from pycbc.io.hdf import SingleDetTriggers

        sdt = SingleDetTriggers.__new__(SingleDetTriggers)
        sdt.ntriggers = 5
        sdt.mask = None

        # 1. Apply mask with list of indices
        sdt.apply_mask([0, 2, 4])
        self.assertEqual(int(sdt.mask_size), 3)

        # 2. Apply scalar index (such as stat.argsort()[0] in page_snglinfo)
        sdt.apply_mask(np.int64(1))
        self.assertEqual(int(sdt.mask_size), 1)
        self.assertTrue(sdt.mask[2])

        # 3. Apply empty list index on fresh object
        sdt2 = SingleDetTriggers.__new__(SingleDetTriggers)
        sdt2.ntriggers = 5
        sdt2.mask = None
        sdt2.apply_mask([])
        self.assertEqual(int(sdt2.mask_size), 0)


suite = unittest.TestSuite()
suite.addTest(unittest.TestLoader().loadTestsFromTestCase(TestIOHDF))


if __name__ == "__main__":
    results = unittest.TextTestRunner(verbosity=2).run(suite)
    simple_exit(results)
