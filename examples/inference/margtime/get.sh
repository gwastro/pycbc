set -e

for ifo in H-H1 L-L1
do
    file=${ifo}_GWOSC_4KHZ_R1-1126257415-4096.gwf
    test -f ${file} && continue
    curl -O --show-error --silent https://www.gw-openscience.org/eventapi/html/GWTC-1-confident/GW150914/v3/${ifo}_GWOSC_4KHZ_R1-1126257415-4096.gwf
done
