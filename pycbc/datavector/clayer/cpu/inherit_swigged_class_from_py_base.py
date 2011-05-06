
import sys,re

filetopatch = sys.stdin.readlines();

cnt = 0
for l in filetopatch:
    l = l.rstrip("\n")
    if re.match("^class %s\(_object\):$" % sys.argv[1], l):
        print "from %s import %s" % (sys.argv[2], sys.argv[3])
        print "class %s(%s, _object):" % (sys.argv[1], sys.argv[3])
        cnt += 1
    else:
        print l

sys.exit(cnt != 1)

