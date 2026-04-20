"""
Test that the constants provided from LAL match closely to the options that
we fall back on from astropy or numpy as appropriate
"""
import unittest
import unittest.mock
import sys
import numpy as np


def get_constant_names(module):
    """Dynamically find uppercase constants, excluding module internals."""
    return [
        c for c in dir(module)
        if c.upper() == c
        and not c.startswith('_') # Exclude all internal vars
    ]

# Define the module name - this is just to make it easier to rename later
_PYCBC_CONSTANTS_MODULE_NAME = 'pycbc.constants'

class TestPycbcConstants(unittest.TestCase):
    """Tests the choice mechanism of pycbc.constants."""

    def setUp(self):
        # We need a fresh import, so we must remove it from cache if it exists.
        if _PYCBC_CONSTANTS_MODULE_NAME in sys.modules:
            del sys.modules[_PYCBC_CONSTANTS_MODULE_NAME]

        # Default constants
        with unittest.mock.patch.dict('os.environ', {'PYCBC_CONSTANT_SOURCE': 'default'}):
            import pycbc.constants as default_constants
        self.default_constants = default_constants

        # Reset the sys.modules value
        if _PYCBC_CONSTANTS_MODULE_NAME in sys.modules:
            del sys.modules[_PYCBC_CONSTANTS_MODULE_NAME]

        # LAL is used here
        with unittest.mock.patch.dict('os.environ', {'PYCBC_CONSTANT_SOURCE': 'lal'}):
            import pycbc.constants as lal_constants
        self.lal_constants = lal_constants

        # Dynamically get the list of constants to check
        self.const_name_list = get_constant_names(default_constants)

    def test_constants_choice(self):
        """Tests that constants are chosen correctly given the environment variable"""
        
        # Assert _CONSTANTS value check for default version
        self.assertTrue(
            self.default_constants._CONSTANTS == 'default',
            "Import has provided wrong constants given environment variable."
        )

        # Assert _CONSTANTS value check for LAL version
        self.assertTrue(
            self.lal_constants._CONSTANTS == 'lal',
            "Import has provided wrong constants given environment variable."
        )

    def test_lal_default_match(self):
        """Tests that LAL constants match default constants."""
        # Use unittest's assertion methods instead of if/print/assert
        for const_name in self.const_name_list:
            lc = getattr(self.lal_constants, const_name)
            fc = getattr(self.default_constants, const_name)

            with self.subTest(constant=const_name):
                # Check for numerical closeness (standard for float comparisons)
                self.assertTrue(
                    np.isclose(lc, fc),
                    msg=(
                        f"Constant {const_name} value mismatch. "
                        f"LAL:{lc}, DEFAULT:{fc}"
                    )
                )

if __name__ == '__main__':
    unittest.main()