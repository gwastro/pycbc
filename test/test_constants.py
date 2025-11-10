"""
Test that the constants provided from LAL match closely to the options that
we fall back on from astropy or numpy as appropriate
"""
import unittest
import sys
from unittest import mock
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
    """Tests the fallback mechanism of pycbc.constants."""



    def test_constants_fallback(self):
        """Tests that constants fall back to Astropy/Numpy when LAL is unavailable."""
        
        # We need a fresh import, so we must remove it from cache if it exists.
        if _PYCBC_CONSTANTS_MODULE_NAME in sys.modules:
            del sys.modules[_PYCBC_CONSTANTS_MODULE_NAME]

        # Default constants
        with mock.patch.dict('os.environ', {'PYCBC_CONSTANT_SOURCE': 'default'}):
            import pycbc.constants as default_constants

        # Assert _CONSTANTS value check for default version
        self.assertTrue(
            default_constants._CONSTANTS == 'default', 
            "Fallback import failed to detect LAL unavailability."
        )

        # Reset the sys.modules value
        if _PYCBC_CONSTANTS_MODULE_NAME in sys.modules:
            del sys.modules[_PYCBC_CONSTANTS_MODULE_NAME]

        # LAL is used here
        with mock.patch.dict('os.environ', {'PYCBC_CONSTANT_SOURCE': 'lal'}):
            import pycbc.constants as lal_constants
        
        # Assert _CONSTANTS value check for LAL version
        self.assertTrue(
            lal_constants._CONSTANTS == 'lal', 
            "LAL-based import failed to detect LAL availability."
        )

        # Dynamically get the list of constants to check
        const_name_list = get_constant_names(default_constants)

        # Use unittest's assertion methods instead of if/print/assert
        for const_name in const_name_list:
            print(f"Checking {const_name}")
            lc = getattr(lal_constants, const_name)
            fc = getattr(default_constants, const_name)

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