set -e

for ifo in H-H1 L-L1 V-V1
do
    file=${ifo}_LOSC_CLN_4_V1-1187007040-2048.gwf
    test -f ${file} && continue

    # Not downloading frames from GWOSC to avoid failure to respond errors on
    # the test suite.
    # The commented command below is how to get the frame from GWOSC if you
    # wanted to verify they are the same.

#    curl -O -L --show-error --silent https://dcc.ligo.org/public/0146/P1700349/001/${file}

    curl -O -L --show-error --silent \
         https://media.githubusercontent.com/media/gwastro/pycbc_data/master/${file}
done
