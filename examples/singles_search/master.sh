#!/bin/bash
set -e

bash -e ../search/get.sh
bash -e ../search/bank.sh
bash -e ../search/stats.sh
bash -e gen.sh

cp *.gwf output
cd output
bash -e ../../search/submit.sh
python ../../search/check_job.py
