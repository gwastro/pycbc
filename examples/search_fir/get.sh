#!/bin/bash
set -e

FRAME_DIR="${PYCBC_FRAME_DIR:-$(pwd)/frames}"
mkdir -p "${FRAME_DIR}"

BASE_URL="https://gwosc.org/archive/data/O2_4KHZ_R1/1185972224"

download_frame() {
    local fname="$1"
    if [ ! -s "${FRAME_DIR}/${fname}" ]; then
        wget -nv "${BASE_URL}/${fname}" -O "${FRAME_DIR}/${fname}"
    fi
}

# Download O2 4 kHz GWOSC frames spanning single-, dual-, and triple-detector
# time (GPS 1186725888 - 1186754560):
#   1186725888-1186729984: 1-det H1
#   1186729984-1186742272: 3-det H1L1V1 + natural 2-det H1L1 & H1V1
#   1186742272-1186746368: 2-det L1V1
#   1186746368-1186750464: 1-det L1
#   1186750464-1186754560: 1-det V1
H1_TIMES=(1186725888 1186729984 1186734080 1186738176)
L1_TIMES=(1186729984 1186734080 1186738176 1186742272 1186746368)
V1_TIMES=(1186729984 1186734080 1186738176 1186742272 1186750464)

: > h1.cache
for t in "${H1_TIMES[@]}"; do
    fn="H-H1_GWOSC_O2_4KHZ_R1-${t}-4096.gwf"
    download_frame "${fn}"
    echo "H H1_GWOSC_O2_4KHZ_R1 ${t} 4096 file://$(realpath "${FRAME_DIR}/${fn}")" >> h1.cache
done

: > l1.cache
for t in "${L1_TIMES[@]}"; do
    fn="L-L1_GWOSC_O2_4KHZ_R1-${t}-4096.gwf"
    download_frame "${fn}"
    echo "L L1_GWOSC_O2_4KHZ_R1 ${t} 4096 file://$(realpath "${FRAME_DIR}/${fn}")" >> l1.cache
done

: > v1.cache
for t in "${V1_TIMES[@]}"; do
    fn="V-V1_GWOSC_O2_4KHZ_R1-${t}-4096.gwf"
    download_frame "${fn}"
    echo "V V1_GWOSC_O2_4KHZ_R1 ${t} 4096 file://$(realpath "${FRAME_DIR}/${fn}")" >> v1.cache
done
