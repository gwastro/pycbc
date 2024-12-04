#!/bin/bash
set -e

if [[ ! -f bank.hdf ]]
then
  bash -e ../../search/bank.sh
else
  echo -e "\\n\\n>> [`date`] Pre-existing template bank found"
fi

bash -e gen.sh

cd output

bash -e ../../../search/submit.sh
python ../../../search/check_job.py

