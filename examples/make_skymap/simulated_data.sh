set -e

# generate an injection file using the PyCBC Live example
../live/generate_injections.py

pycbc_make_skymap \
    --trig-time 1272790260 \
    --fake-strain \
        H1:aLIGOMidLowSensitivityP1200087 \
        L1:aLIGOMidLowSensitivityP1200087 \
        V1:AdVEarlyHighSensitivityP1200087 \
    --injection-file injections.hdf \
    --thresh-SNR 5.5 \
    --f-low 20 \
    --mass1 1.1331687 \
    --mass2 1.010624 \
    --spin1z 0.029544285 \
    --spin2z 0.020993788 \
    --ifos H1 L1 V1 \
    --ligolw-event-output coinc_simulated_data.xml

# make a zoomed-in plot to see the small skymap and the true position
ligo-skymap-plot \
    --output 1272790260_skymap_zoom.png \
    --projection zoom \
    --projection-center '10deg -45deg' \
    --zoom-radius 20deg \
    --contour 50 90 \
    --annotate \
    --radec 10 -45 \
    1272790260.fits
