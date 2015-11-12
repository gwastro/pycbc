#!/bin/bash

pushd ../../pycbc/workflow 

PEGASUS_FILE_DATA=`tar -zc pegasus_files | base64 | xargs echo | sed 's+ +ZZZZ+g'`

popd 

mkdir -p dist

sed -e 's+DATA_INLINE=False+DATA_INLINE=True+' \
    -e "s^PEGASUS_FILE_DATA^${PEGASUS_FILE_DATA}^" \
    -e "s+ZZZZ+\n+g" ../../bin/pycbc_submit_dax > dist/pycbc_submit_dax


