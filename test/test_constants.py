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
        and not c.startswith('_') # Exclude all dunders and internal vars
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

        # LAL is available here (imported before the test function or via setup, 
        # but Python finds it as a true import)
        import pycbc.constants as lal_constants
        
        # Assert availability check is correct for the LAL version
        # i.e. if we run this test without lal being available,
        # this check will fail
        self.assertTrue(lal_constants._LAL_AVAILABLE, 
                        "LAL-based import failed to detect LAL availability.")

        # Now load to get the fallback constants
        # Re-remove from cache *before* the second import.
        # This is critical because the @mock.patch only affects 'lal', 
        # not 'pycbc.constants'.
        if _PYCBC_CONSTANTS_MODULE_NAME in sys.modules:
            del sys.modules[_PYCBC_CONSTANTS_MODULE_NAME]

        with mock.patch.dict('sys.modules', {'lal': None}):
            import pycbc.constants as fallback_constants

        # Assert availability check is correct for the fallback version
        self.assertFalse(fallback_constants._LAL_AVAILABLE, 
                         "Fallback import failed to detect LAL unavailability.")

        # --- Stage 3: Compare Constants ---
        
        # Dynamically get the list of constants to check
        const_name_list = get_constant_names(lal_constants)

        # Use unittest's assertion methods instead of if/print/assert
        for const_name in const_name_list:
            lc = getattr(lal_constants, const_name)
            fc = getattr(fallback_constants, const_name)

            with self.subTest(constant=const_name):
                # Check for numerical closeness (standard for float comparisons)
                self.assertTrue(
                    np.isclose(lc, fc),
                    msg=(
                        f"Constant {const_name} value mismatch. "
                        f"LAL:{lc}, FALLBACK:{fc}"
                    )
                )
                
# Optional: Boilerplate to run the test if the file is executed directly
if __name__ == '__main__':
    unittest.main()