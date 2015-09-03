######################
pycbc_inspiral
######################

Overview
========

This code takes in gravitational wave data, matches it against a bank of known waveforms, and returns a list of triggers in a LIGO_LW format XML file.

Help message
============

.. command-output:: pycbc_inspiral --help 

GPU
===

Instructions on how to get pycbc_inspiral working on GPUs. Currently, these instructions only work for CUDA 6.5 installation. You will need to install a particular version of pycbc. Currently, this version works:

.. code-block:: bash

    pip install git+https://github.com/ligo-cbc/pycbc@f5f1775a5d6587066d0d1ad0b9e27a52e8aef68c  --process-dependency-links

Check to make sure that you have pyfft.cuda installed. If you don't, you will also need to install it:

.. code-block:: bash

    pip install https://pypi.python.org/packages/source/p/pyfft/pyfft-0.3.tar.gz

Hopefully, you have installed everything and it was relatively painless (if it wasn't, you have my deepest sympathies). But let's move on!

Activate your virtual environment, and create a directory for testing. 

.. code-block:: bash
    
    . activate.sh # Script that activates your virtual environment
    mkdir test

Inside this directory, you will set up three different data sets (clean, loud, grumbly).

.. code-block:: bash

    cd test
    
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

Copy over the template bank for testing. Start with the small template bank for initial testing, and if that works, then try the large template bank.

.. code-block:: bash

    cp /home/lppekows/Profiling/TMPLTBANK_SMALL.xml.gz . # Small template bank
    cp /home/lppekows/Profiling/TMPLTBANK.xml.gz . # Large template bank

Now, run the actual pycbc_inspiral command and its several inputs:

.. code-block:: bash

    pycbc_inspiral 
    --cluster-method window \
    --cluster-window 1  \
    --bank-file TMPLTBANK.xml.gz \
    --approximant SPAtmplt  \
    --gps-start-time 970910372  \
    --gps-end-time   970912420  \
    --snr-threshold 5.5  \
    --strain-high-pass 25.0  \
    --chisq-bins 16  \
    --psd-inverse-length 16  \
    --psd-segment-stride 128  \
    --psd-segment-length 256  \
    --psd-estimation median  \
    --segment-length 256  \
    --segment-start-pad 112  \
    --segment-end-pad 16  \
    --low-frequency-cutoff 30.0  \
    --pad-data 8  \
    --sample-rate 4096  \
    --order 7  \
    --frame-files clean/*.gwf  \
    --channel-name H1:LDAS-STRAIN  \
    --output test_0-5.hdf5  \
    --processing-scheme cuda

You will need to specify the following, as they may vary between different tests/data sets (others can be left as default):

.. code-block:: bash

    --bank-file # The template bank you will use (specify directory if needed)
    --gps-start-time # Start time of run
    --gps-end-time # End time of run
    --frame-files # The directory where the frame files are
    --output # Name of the output file (specify directory if needed)
    
Here is the GPS_START_TIME and GPS_END_TIME for each data set:

* clean: GPS_START_TIME: 970910372, GPS_END_TIME: 970912420
* loud: GPS_START_TIME: 962784401, GPS_END_TIME: 962786449
* grumbly: GPS_START_TIME: 969658210, GPS_END_TIME: 969660258

If you have virtualenv and all other software installed properly, you should be able to run the pycbc_inspiral command and have it go smoothly. 

Optional: If you choose to use screen (so that you can logout and log back in later/on a different computer), run the commands in the following order:

.. code-block:: bash
  
    screen
    . activate.sh # Script that activates your virtual environment
    cd path/to/frame-files/ # Go to where your frame-files directory is
    pycbc_inspiral ... # pycbc_inspiral, with all necessary inputs
    
Running in a different order may cause an error (example: "The 'PyCBC===4feb06' distribution was not found and is required by the application").
