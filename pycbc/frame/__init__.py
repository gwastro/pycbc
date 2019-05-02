from . frame import (locations_to_cache, read_frame, datafind_connection,
                     query_and_read_frame, frame_paths, write_frame,
                     DataBuffer, StatusBuffer)

from . store import (read_store)


# Status flags for the calibration state vector
# See e.g. https://dcc.ligo.org/LIGO-G1700234
# https://wiki.ligo.org/DetChar/DataQuality/O3Flags
HOFT_OK = 1
SCIENCE_INTENT = 2
SCIENCE_QUALITY = 4
HOFT_PROD = 8
FILTERS_OK = 16
NO_STOCH_HW_INJ = 32
NO_CBC_HW_INJ = 64
NO_BURST_HW_INJ = 128
NO_DETCHAR_HW_INJ = 256
KAPPA_A_OK = 512
KAPPA_PU_OK = 1024
KAPPA_TST_OK = 2048
KAPPA_C_OK = 4096
FCC_OK = 8192
NO_GAP = 16384
NO_HWINJ = NO_STOCH_HW_INJ | NO_CBC_HW_INJ | \
           NO_BURST_HW_INJ | NO_DETCHAR_HW_INJ

# relevant bits in the LIGO O2/O3 low-latency DQ vector
# If the bit is 0 then we should veto
# https://wiki.ligo.org/DetChar/DmtDqVector
# https://wiki.ligo.org/DetChar/DataQuality/O3Flags
OMC_DCPD_ADC_OVERFLOW = 2
ETMY_ESD_DAC_OVERFLOW = 4
ETMX_ESD_DAC_OVERFLOW = 16

# CAT1 bit in the Virgo state vector
# https://wiki.virgo-gw.eu/DetChar/DetCharVirgoStateVector
VIRGO_GOOD_DQ = 1 << 10


def flag_names_to_bitmask(flags):
    """Takes a list of flag names corresponding to bits in a status channel
    and returns the corresponding bit mask.
    """
    mask = 0
    for flag in flags:
        mask |= globals()[flag]
    return mask
