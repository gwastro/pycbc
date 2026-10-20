#!/bin/bash
# Install HTCondor and Pegasus on Ubuntu CI runners; requires sudo.
set -euo pipefail

# Configure the HTCondor repository and install the local scheduler.
wget -qO - https://research.cs.wisc.edu/htcondor/ubuntu/HTCondor-Release.gpg.key | sudo apt-key add -
echo "deb http://research.cs.wisc.edu/htcondor/ubuntu/8.9/focal focal contrib" | sudo tee -a /etc/apt/sources.list
echo "deb-src http://research.cs.wisc.edu/htcondor/ubuntu/8.9/focal focal contrib" | sudo tee -a /etc/apt/sources.list
sudo apt-get -o Acquire::Retries=3 update
sudo apt-get -o Acquire::Retries=3 install -y minihtcondor

# Configure the Pegasus repository and install the workflow manager.
wget -qO - https://download.pegasus.isi.edu/pegasus/gpg.txt | sudo apt-key add -
echo "deb https://download.pegasus.isi.edu/pegasus/ubuntu noble main" | sudo tee -a /etc/apt/sources.list
sudo apt-get -o Acquire::Retries=3 update
sudo apt-get -o Acquire::Retries=3 install -y pegasus=5.1.1-1+ubuntu24
