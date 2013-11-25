#set path to here
PATH=$PWD:$PATH

#get frames
sh df.sh

#run pycbc_inspiral
sh prun.sh

#run lalapps_inspira
sh lrun.sh

# run trigger comparison
python trig-compare.py H1-INSPIRAL_pycbc_FULL_DATA-968605000-2048.xml.gz H1-INSPIRAL_lalsuite_FULL_DATA-968605000-2048.xml.gz
