#!/bin/bash
set -e

GWDATA_DIR="${PYCBC_GWDATA_DIR:-}"
FRAME_DIR="${PYCBC_FRAME_DIR:-$(pwd)/frames}"
BASE_URL="https://gwosc.org/archive/data/O2_4KHZ_R1/1185972224"

download_frame() {
    local fname="$1"
    mkdir -p "${FRAME_DIR}"
    if [ ! -s "${FRAME_DIR}/${fname}" ]; then
        wget -nv "${BASE_URL}/${fname}" -O "${FRAME_DIR}/${fname}"
    fi
}

# Build HDF5 frame caches spanning single-, dual-, and triple-detector
# time (GPS 1186725888 - 1186754560):
#   1186725888-1186729984: 1-det H1
#   1186729984-1186742272: 3-det H1L1V1 + natural 2-det H1L1 & H1V1
#   1186742272-1186746368: 2-det L1V1
#   1186746368-1186750464: 1-det L1
#   1186750464-1186754560: 1-det V1
H1_TIMES=(1186725888 1186729984 1186734080 1186738176)
L1_TIMES=(1186729984 1186734080 1186738176 1186742272 1186746368)
V1_TIMES=(1186729984 1186734080 1186738176 1186742272 1186750464)

write_cache() {
    local site="$1"
    local ifo="$2"
    local out="$3"
    shift 3
    : > "${out}"
    for t in "$@"; do
        local local_2k="${GWDATA_DIR}/O2/${ifo}/1185939456/${site}-${ifo}_DATASAMPLER_O2_2048HZ-${t}-4096.hdf5"
        if [ -s "${local_2k}" ]; then
            echo "${site} ${ifo}_DATASAMPLER_O2_2048HZ ${t} 4096 file://$(realpath "${local_2k}")" >> "${out}"
        else
            local fn="${site}-${ifo}_GWOSC_O2_4KHZ_R1-${t}-4096.hdf5"
            download_frame "${fn}"
            echo "${site} ${ifo}_GWOSC_O2_4KHZ_R1 ${t} 4096 file://$(realpath "${FRAME_DIR}/${fn}")" >> "${out}"
        fi
    done
}

write_cache H H1 h1.cache "${H1_TIMES[@]}"
write_cache L L1 l1.cache "${L1_TIMES[@]}"
write_cache V V1 v1.cache "${V1_TIMES[@]}"
