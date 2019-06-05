# Custom environment shell code should follow

# If no LAL_DATA_PATH has been set, use the data in the container
if [ "x$LAL_DATA_PATH" == "x" ]; then
    export LAL_DATA_PATH="/cvmfs/oasis.opensciencegrid.org/ligo/sw/pycbc/lalsuite-extra/current/share/lalsimulation"
fi

# Add path to Intel MKL Libraries
if [ "x$LD_LIBRARY_PATH" == "x" ]; then
    LD_LIBRARY_PATH="/opt/intel/composer_xe_2015.0.090/mkl/lib/intel64:/opt/mvapich2-2.1/lib"
else
    LD_LIBRARY_PATH="/opt/intel/composer_xe_2015.0.090/mkl/lib/intel64:/opt/mvapich2-2.1/lib:${LD_LIBRARY_PATH}"
fi
export LD_LIBRARY_PATH

PATH=${PATH}:/opt/mvapich2-2.1/bin
export PATH
