set -e

pycbc_make_skymap \
    --trig-time 1187008882.4457 \
    --thresh-SNR 5.5 \
    --f-low 27 \
    --mass1 1.457423 \
    --mass2 1.299302 \
    --spin1z "-0.018811" \
    --spin2z "0.011777" \
    --ifos H1 L1 V1 \
    --ligolw-event-output coinc_GW170817.xml

# make a zoomed-in plot to see the details
ligo-skymap-plot \
    --output 1187008882_skymap_zoom.png \
    --projection zoom \
    --projection-center '13h8m -21d' \
    --zoom-radius 10deg \
    --contour 50 90 \
    --annotate \
    --radec 197.4500 -23.3814 \
    1187008882.fits
