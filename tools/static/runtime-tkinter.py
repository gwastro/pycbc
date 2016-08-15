import os, sys
d = os.path.join(sys._MEIPASS, 'tcl')
if not os.path.exists(d):
    os.makedirs(d)
d = os.path.join(sys._MEIPASS, 'tk')
if not os.path.exists(d):
    os.makedirs(d)
