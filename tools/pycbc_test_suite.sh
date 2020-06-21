#!/bin/bash

echo -e "\\n>> [`date`] Starting PyCBC test suite"

LOG_FILE=$(mktemp -t pycbc-test-log.XXXXXXXXXX)

RESULT=0

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
