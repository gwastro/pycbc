#!/bin/bash
set -e

echo -e ">> [`date`] upgrading setuptools and pip"
pip install --upgrade setuptools pip

echo -e ">> [`date`] installing requirements"
pip install -r requirements.txt
