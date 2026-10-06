#!/bin/bash
# Start HTCondor and wait until its scheduler (schedd) is reachable.
#
# minihtcondor must already be installed.  Used by the pegasus + condor CI
# workflows.  On a fresh runner condor.service sometimes comes up "active" but
# with no daemons actually running (condor_master is alive, but there is no
# schedd); a workflow submit then fails later with the cryptic
#   ERROR: Can't find address of local schedd
# so here we make the hostname resolvable, start condor, and restart it a few
# times until `condor_q` succeeds, dumping the condor logs if it never does.

# The condor daemons fail to start if the machine's own hostname does not
# resolve.
if ! grep -q "$(hostname)" /etc/hosts; then
    echo "127.0.1.1 $(hostname)" | sudo tee -a /etc/hosts
fi

sudo systemctl enable condor

for start in 1 2 3; do
    sudo systemctl restart condor
    for poll in $(seq 1 20); do
        if condor_q > /dev/null 2>&1; then
            echo "condor schedd ready (start ${start}, after ${poll} polls)"
            condor_status
            exit 0
        fi
        echo "waiting for condor schedd (start ${start}, poll ${poll}/20)..."
        sleep 3
    done
    echo "condor schedd still not up after start ${start}"
    sudo systemctl status condor --no-pager
    sudo tail -n 50 /var/log/condor/MasterLog 2>/dev/null
done

echo "::error::condor schedd did not become ready"
sudo systemctl status condor --no-pager
sudo journalctl -u condor --no-pager | tail -n 50
sudo tail -n 100 /var/log/condor/MasterLog /var/log/condor/SchedLog 2>/dev/null
exit 1
