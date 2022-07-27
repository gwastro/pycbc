# About 482 templates

pycbc_geom_nonspinbank \
    --pn-order twoPN \
    --f0 40 \
    --f-low 40 \
    --f-upper 2048 \
    --delta-f 0.1 \
    --min-match 0.97 \
    --min-mass1 2.0 \
    --min-mass2 2.0 \
    --max-mass1 3 \
    --max-mass2 3 \
    --verbose \
    --output-file testNonSpin2.xml \
    --calculate-ethinca-metric \
    --filter-cutoff SchwarzISCO \
    --asd-file ZERO_DET_high_P.txt

ligolw_print \
    -t sngl_inspiral \
    -c mass1 \
    -c mass2 \
    testNonSpin2.xml
