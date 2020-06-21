#!/bin/bash

echo -e "\\n>> [`date`] Starting PyCBC test suite"

LOG_FILE=$(mktemp -t pycbc-test-log.XXXXXXXXXX)

RESULT=0

echo "MYTEST"

python -c "import numpy
a = numpy.array([-2.8498888-4.8523064j,
                          11.672292 -0.39886224j,
                           -3.907175 -3.9876463j,
                           3.5721998-4.4356446j,
                           -0.9019761+6.379277j], dtype=numpy.complex64)
print (abs(a)) 
"

#run pycbc inspiral test
pushd examples/inspiral
bash -e run.sh
if test $? -ne 0 ; then
    RESULT=1
    echo -e "    FAILED!"
    echo -e "---------------------------------------------------------"
else
    echo -e "    Pass."
fi
popd

exit ${RESULT}
