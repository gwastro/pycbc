"""Check template parameters and frequency metadata survive bank output."""
import tempfile
import unittest
from pathlib import Path

import h5py
import numpy
from igwn_ligolw import lsctables, utils

from pycbc.io.ligolw import LIGOLWContentHandler
from pycbc.tmpltbank.bank_output_utils import output_bank_to_file


class BankOutputTest(unittest.TestCase):
    """Exercise both bank formats without requiring PSD downloads."""

    def test_output_roundtrip(self):
        templates = [(3.0, 2.0, 0.1, -0.05), (4.0, 1.5, -0.2, 0.0)]
        options = {"f_low": 30.0, "f_upper": 512.0, "output_f_final": True,
                   "f_low_column": 'alpha'}
        with tempfile.TemporaryDirectory() as directory:
            for suffix in ['xml', 'xml.gz', 'hdf']:
                with self.subTest(suffix=suffix):
                    path = str(Path(directory) / ('bank.' + suffix))
                    output_bank_to_file(path, templates, opt_dict=options)
                    if suffix == 'hdf':
                        with h5py.File(path, 'r') as bank:
                            for i, column in enumerate(
                                    ['mass1', 'mass2', 'spin1z', 'spin2z']):
                                numpy.testing.assert_array_equal(
                                    bank[column][:],
                                    numpy.array(templates)[:, i])
                            numpy.testing.assert_array_equal(
                                bank['f_lower'][:], [30.0, 30.0])
                            numpy.testing.assert_array_equal(
                                bank['f_final'][:], [512.0, 512.0])
                    else:
                        document = utils.load_filename(
                            path, contenthandler=LIGOLWContentHandler)
                        bank = lsctables.SnglInspiralTable.get_table(document)
                        self.assertEqual(len(bank), len(templates))
                        for row, expected in zip(bank, templates, strict=True):
                            numpy.testing.assert_allclose(
                                [row.mass1, row.mass2, row.spin1z, row.spin2z],
                                expected, rtol=1e-6)
                            self.assertEqual(row.alpha, 30.0)
                        document.unlink()


if __name__ == '__main__':
    unittest.main()
