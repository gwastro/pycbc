#!/bin/bash
set -e

bash -e get.sh
bash -e bank.sh
bash -e stats.sh
bash -e generate_unknown_injections.sh
bash -e gen.sh

cp *.gwf output
cd output
bash -e ../submit.sh
python ../check_job.py
