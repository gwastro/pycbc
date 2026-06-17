set -e

for ifo in H-H1 L-L1
do
    file=${ifo}_GWOSC_4KHZ_R1-1126257415-4096.gwf
    test -f ${file} && continue

    # Not downloading frames from GWOSC to avoid failures.
    # GWOSC often is not responsive to queries from within the GitHub CI.
    # The commented command below is how to get the frame from GWOSC if you
    # wanted to verify they are the same.

    #curl -O -L --show-error --silent \
    #    https://www.gwosc.org/eventapi/html/GWTC-1-confident/GW150914/v3/${file}
    curl -O -L --show-error --silent \
	 https://media.githubusercontent.com/media/gwastro/pycbc_data/master/${file}
done
