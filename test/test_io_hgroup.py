import os
import unittest
import tempfile
import h5py
import numpy as np
from utils import simple_exit, parse_args_cpu_only
from pycbc.io.hdf import HFile, HGroup

parse_args_cpu_only('io.hgroup')


def _dataset_has_fletcher32(dataset):
    """Return True if dataset creation property lists include fletcher32."""
    plist = dataset.id.get_create_plist()
    n = plist.get_nfilters()
    for i in range(n):
        tpl = plist.get_filter(i)
        name = tpl[-1]
        if isinstance(name, bytes):
            name_low = name.lower()
        else:
            name_low = str(name).lower()
        if b'fletcher' in name_low if isinstance(name_low, bytes) else 'fletcher' in name_low:
            return True
    return False


class TestIOHFile(unittest.TestCase):
    def test_create_dataset_has_fletcher32(self):
        """
        Test that a created dataset has fletcher32 enabled or not as appropriate"""
        should_be_checksummed = []
        should_not_be_checksummed = []
        with tempfile.TemporaryDirectory() as td:
            # Lets make a file / some datasets using HFile:
            p = os.path.join(td, "f1.hdf")
            with HFile(p, 'w') as f:
                # This is a dataset which has been initialised without any data,
                # but has a numeric dtype set
                f.create_dataset(
                    'shape_based_with_dtype',
                    shape=(10,),
                    dtype='f8'
                )
                should_be_checksummed.append('shape_based_with_dtype')

                # Make a shape-based dataset where dtype isnt given
                f.create_dataset(
                    'shape_based_no_dtype',
                    shape=(10,),
                )
                should_not_be_checksummed.append('shape_based_no_dtype')

                int_list = [1,2,3,4,5,6,7,8,9]
                # Make a dataset from some int data cast as a f8,
                f.create_dataset(
                    'data_np_int_f8_array',
                    data=np.array(int_list, dtype=np.int32),
                    dtype='f8'
                )
                should_be_checksummed.append('data_np_int_f8_array')

                # Make a dataset from some int data,
                f.create_dataset(
                    'data_np_int_array',
                    data=np.array(
                        int_list,
                        dtype=np.int32
                    ),
                )
                should_be_checksummed.append('data_np_int_array')

                # Make a dataset from an int list,
                f.create_dataset(
                    'data_int_list',
                    data=int_list,
                )
                should_be_checksummed.append('data_int_list')

                str_list = ['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h']
                # Make a dataset from a string numpy array data
                f.create_dataset(
                    'data_np_str_array',
                    data=np.array(str_list, dtype='S1')
                )
                should_not_be_checksummed.append('data_np_str_array')

                # Make a dataset from some string data
                f.create_dataset(
                    'data_str_list',
                    data=str_list,
                    dtype=h5py.string_dtype('utf-8')
                )
                should_not_be_checksummed.append('data_str_list')

                # Make an object array from some string data
                f.create_dataset(
                    'data_obj_list',
                    data=np.array(str_list, dtype=object)
                )
                should_not_be_checksummed.append('data_obj_list')

                # Make a datasets from scalars:
                f.create_dataset(
                    'data_scalar_with_dtype',
                    data=0,
                    dtype='f8'
                )
                should_not_be_checksummed.append('data_scalar_with_dtype')
                f.create_dataset(
                    'data_scalar_no_dtype',
                    data=0,
                )
                should_not_be_checksummed.append('data_scalar_no_dtype')

            with h5py.File(p, 'r') as f:
                for dataset_name in should_be_checksummed:
                    d = f[dataset_name]
                    self.assertTrue(
                        _dataset_has_fletcher32(d),
                        "fletcher32 should be enabled for %s dataset" % dataset_name
                    )
                for dataset_name in should_not_be_checksummed:
                    d = f[dataset_name]
                    self.assertFalse(
                        _dataset_has_fletcher32(d),
                        "fletcher32 should be not enabled for %s dataset" % dataset_name
                    )

    def test_getitem_returns_hgroup_and_root_parent_behaviour(self):
        with tempfile.TemporaryDirectory() as td:
            p = os.path.join(td, "f3.hdf")
            with HFile(p, 'w') as f:
                grp = f.create_group('grp')
                self.assertIsInstance(grp, HGroup,
                                      f"Expected group to be HGroup, got {type(grp)}")
                subgrp = grp.create_group('subgrp')
                self.assertIsInstance(subgrp, HGroup,
                                      f"Expected subgroup to be HGroup, got {type(subgrp)}")


                root = f['/']
                self.assertIsInstance(root.parent, HGroup,
                                      "Expected root.parent to be an HGroup")
                self.assertEqual(root.parent.id, root.id,
                                 "Expected root.parent to refer to the same group as root")

    def test_create_dataset_incompatible_dtype_raises(self):
        """Creating a dataset where data cannot be converted to the
        requested dtype should raise a TypeError or ValueError.
        """
        with tempfile.TemporaryDirectory() as td:
            p = os.path.join(td, "bad_dtype.hdf")
            # create a unicode numpy array whose elements are non-numeric
            arr = np.array(['a', 'b', 'c'], dtype='<U1')
            with HFile(p, 'w') as f:
                # Attempting to write this array as floats should fail
                with self.assertRaises((TypeError, ValueError)):
                    f.create_dataset('bad', data=arr, dtype='f8')

    def test_getitem_wraps_group_and_leaves_dataset(self):
        """
        Ensure __getitem__ returns HGroup for groups and leaves datasets
        untouched as datasets.
        """
        with tempfile.TemporaryDirectory() as td:
            p = os.path.join(td, "f5.hdf")
            with HFile(p, 'w') as f:
                f.create_group('grp2')
                # create a simple dataset
                f.create_dataset('ds_plain', data=np.arange(3, dtype=np.int32))

            # Open with plain h5py to test what is returned through normal access
            with HFile(p, 'r') as f:
                got_grp = f['grp2']
                got_ds = f['ds_plain']


                self.assertIsInstance(
                    got_grp,
                    HGroup,
                    f"Expected a HGroup, got {type(got_grp)}"
                )

                self.assertIsInstance(
                    got_ds,
                    h5py.Dataset,
                    f"Expected a h5py.Dataset, got {type(got_ds)}"
                )


    def test_dataset_parity_h5py_vs_hfile(self):
        """Create identical datasets with h5py.File and HFile and compare parity."""
        with tempfile.TemporaryDirectory() as td:
            p_hfile = os.path.join(td, "hf.hdf")
            p_h5py = os.path.join(td, "hp.hdf")

            data = np.arange(20, dtype=np.float64)

            # create with HFile (our wrapper)
            with HFile(p_hfile, 'w') as f:
                f.create_dataset('ds', data=data, dtype='f8')

            # create with plain h5py
            with h5py.File(p_h5py, 'w') as f:
                f.create_dataset('ds', data=data, dtype='f8')

            # reopen with h5py to inspect the raw datasets
            with h5py.File(p_hfile, 'r') as fh, h5py.File(p_h5py, 'r') as fp:
                d_hfile = fh['ds']
                d_h5py = fp['ds']

                # Data parity
                np.testing.assert_array_equal(
                    d_hfile[:],
                    d_h5py[:],
                    err_msg="Dataset contents should be equal"
                )

                # dtype and shape parity
                self.assertEqual(
                    d_hfile.shape,
                    d_h5py.shape,
                    "Shapes should be identical"
                )
                self.assertEqual(
                    str(d_hfile.dtype),
                    str(d_h5py.dtype),
                    "Dtypes should be identical"
                )

                # Assert HFile-created dataset has fletcher32 enabled
                self.assertTrue(
                    _dataset_has_fletcher32(d_hfile),
                    "HFile-created dataset should have fletcher32 checksum"
                )

                # Assert that the h5py.File-created dataset _doesnt_
                # have fletcher32 enabled (If this fails, h5py has updated its
                # default behaviour, and our wrapper isn't needed)
                self.assertFalse(
                    _dataset_has_fletcher32(d_h5py),
                    "h5py.File-created dataset should not have fletcher32 "
                    "checksum (h5py default behaviour has changed - we can remove "
                    "the HGroup wrapper now)"
                )

suite = unittest.TestSuite()
suite.addTest(unittest.TestLoader().loadTestsFromTestCase(TestIOHFile))


if __name__ == '__main__':
    results = unittest.TextTestRunner(verbosity=2).run(suite)
    simple_exit(results)
