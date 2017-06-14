from . frame import (locations_to_cache, read_frame, datafind_connection,
                     query_and_read_frame, frame_paths, write_frame,
                     DataBuffer, StatusBuffer)
                     

# Status flags for the calibration state vector
# See e.g. https://dcc.ligo.org/LIGO-G1700234
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

# O2 Low-Latency DQ vector definition
# If the bit is 0 then we should veto
# https://wiki.ligo.org/DetChar/DmtDqVector
OMC_DCPD_ADC_OVERFLOW = 2
ETMY_ESD_DAC_OVERFLOW = 4

# Virgo state vector
# https://wiki.virgo-gw.eu/DetChar/DetCharVirgoStateVector
VIRGO_GOOD_DQ = 1 << 10
