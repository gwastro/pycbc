#!/usr/bin/env python

import sys

try:
    tottime = float(sys.argv[1])
except:
    print(f"usage: {sys.argv[0]} total_time_in_seconds")
    print()
    print("Typical use case is: ")
    print(
        f"  gprof2dot.py -f pstats profile_file | {sys.argv[0]} total_time_in_seconds "
        "| dot -Tpng -o output.png"
    )

    sys.exit(0)

for line in sys.stdin:
    newtokens = []

    for token in line.strip().split("\\n"):
        if token[-1] == "%":
            try:
                value = float(token[:-1]) / 100.0
                token = f"{token} [{value * tottime:.2e}]"
            except:
                pass
        elif token[-2:] == "%)":
            value = float(token[1:-2]) / 100.0
            token = f"({token[1:-2]} [{value * tottime:.2e}])"

        newtokens.append(token)

    print("\\n".join(newtokens))
