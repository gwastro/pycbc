#!/bin/bash
set -e

bash -e get.sh
bash -e bank.sh
bash -e stats.sh
bash -e generate_unknown_injections.sh
bash -e gen.sh

cd output
ln -sf ../*.gwf .
bash -e ../submit.sh
#python ../check_job.py
