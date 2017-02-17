#!/usr/bin/env python

from __future__ import print_function
import sys
import h5py

datasets = ['template_hash', 'chisq', 'chisq_dof', 'coa_phase', 'end_time', 'sigmasq', 'snr', 'template_duration']

# open file
fin = h5py.File(sys.argv[1],'r')
# print fin.keys()[0]

'''
for group in fin.values():
    print group.keys()
    for dataset in group.values():
        print dataset.name, dataset.len()
print "groups printed"
print fin.items()['/H1/template_hash']
#gin = fin.group()
'''

# should be "/H1/", but also work for e.g. "/L1/"
group = "/" + fin.keys()[0] + "/"

# check that all datasets exist and have same length
l = None
for dataset in datasets:
    if not fin.get(group + dataset):
        raise ValueError("dataset " + group + dataset + " not found")
    elif not l:
        l = fin.get(group + dataset).len()
    elif l != fin.get(group + dataset).len():
        raise ValueError("datasets have different length")

print("Dataset Length:", l, file=sys.stderr)

# print data to stdout in correct order
for i in range(l):
    for dataset in datasets:
        print (" ", fin.get(group + dataset)[i], end="")
    print ("")

fin.close()

sys.exit(0)
