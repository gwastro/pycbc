#!/bin/bash

mkdir clean
cd clean
for i in `ligo_data_find -o H -s 970910336 -e 970923236  --type  H1_LDAS_C02_L2 | grep file | grep -v archive | sed 's+file://localhost++'`
do
	ln -s $i
done
cd ..

mkdir loud
cd loud
for i in `ligo_data_find -o H -s 962784384 -e 962797256  --type  H1_LDAS_C02_L2 | grep file | grep -v archive | sed 's+file://localhost++'`
do
	ln -s $i
done
cd ..

mkdir grumbly
cd grumbly
for i in `ligo_data_find -o H -s 969658112 -e 969670984  --type  H1_LDAS_C02_L2 | grep file | grep -v archive | sed 's+file://localhost++'`
do
	ln -s $i
done
cd ..

