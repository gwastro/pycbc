import subprocess
import time
time.sleep(30)
while 1:
    time.sleep(60)
    subprocess.run(["condor_q", "-bet"])
    subprocess.run(["pegasus-status", "submitdir/work/"])
    out = subprocess.check_output(["pegasus-status", "submitdir/work/"])
    out = str(out)
    lines = out.split('\\n')
    for i in range(len(lines)):
        if 'UNREADY' in lines[i]:
            status_line = i + 1
            break

    stats = lines[status_line].split(' ')
    stats = [s for s in stats if s != '']

    unready = int(stats[0])
    ready = int(stats[1])
    pre = int(stats[2])
    queued = int(stats[3])
    post = int(stats[4])
    done = int(stats[5])
    failed = int(stats[6])

    finished = (unready == 0 and ready == 0 and queued == 0 and post == 0)
    passed = finished and failed == 0

    if passed:
        print("workflow has completed successfully")
        time.sleep(30)
        exit(0)

    if failed != 0:
        print("workflow has a failed job, ending now")
        subprocess.run(["bash", "./stop"])
        # Need to wait here to make sure it fully exits before uploading logs!
        time.sleep(30)
        exit(1)
