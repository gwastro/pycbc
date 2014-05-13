# Non-spinning stochastic bank with ethinca calc
# About 52,000 templates expected
pycbc_aligned_stoch_bank --verbose --pn-order threePointFivePN --f0 15 --f-low 15 --delta-f 0.1 --min-match 0.97 --min-mass1 2. --max-mass1 100. --min-mass2 2. --max-mass2 100. --max-ns-spin-mag 0.0 --max-bh-spin-mag 0.0 --nsbh-flag --verbose --asd-file ZERO_DET_high_P.txt --num-seeds 1000000 --output-file "testStoch_nonspin.xml" --max-total-mass 100 --min-total-mass 10 --max-eta 0.25 --min-eta 0.08264 --calculate-ethinca-metric --f-upper 2000
