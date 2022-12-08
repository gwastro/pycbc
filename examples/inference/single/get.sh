set -e

for ifo in H-H1 L-L1 V-V1
do
    file=${ifo}_LOSC_CLN_4_V1-1187007040-2048.gwf
    test -f ${file} && continue
    curl -O --show-error --silent https://dcc.ligo.org/public/0146/P1700349/001/${file}
done
