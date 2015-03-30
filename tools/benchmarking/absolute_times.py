#!/usr/bin/env python

import sys

try:
    tottime = float(sys.argv[1])
except:
    print "usage: %s total_time_in_seconds" % sys.argv[0]
    print
    print "Typical use case is: "
    print "  gprof2dot.py -f pstats profile_file | %s total_time_in_seconds | dot -Tpng -o output.png"  % sys.argv[0]

    sys.exit(0)

for line in sys.stdin:
    newtokens = []

    for token in line.strip().split('\\n'):
        if token[-1] == '%':
            try:
                value = float(token[:-1]) / 100.0
                token = '%s [%.2e]' % (token, value*tottime)
            except:
                pass
        elif token[-2:] == '%)':
            value = float(token[1:-2]) / 100.0
            token = '(%s [%.2e])' % (token[1:-2], value*tottime)

        newtokens.append(token)

    print '\\n'.join(newtokens)

