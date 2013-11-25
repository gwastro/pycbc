import pylab
from glue.ligolw import table, utils, lsctables
from sys import argv
import numpy

pycbc_file = argv[1]
lal_file = argv[2]

pdoc = utils.load_filename(pycbc_file)
ldoc = utils.load_filename(lal_file)

ptrig = table.get_table(pdoc, lsctables.SnglInspiralTable.tableName)
ltrig = table.get_table(ldoc, lsctables.SnglInspiralTable.tableName)

t = table.get_table(pdoc, lsctables.FilterTable.tableName)
start = t[0].start_time

psnr = ptrig.get_column('snr')
pchi = ptrig.get_column('chisq')
pmc = ptrig.get_column('mchirp')
psig = ptrig.get_column('sigmasq')
pend_time = (ptrig.get_column('end_time') - start)+ 1e-9 * ptrig.get_column('end_time_ns') 

lsnr = ltrig.get_column('snr')
lchi = ltrig.get_column('chisq')
lmc = ltrig.get_column('mchirp')
lsig = ltrig.get_column('sigmasq')
lend_time = (ltrig.get_column('end_time') - start) + 1e-9 * ltrig.get_column('end_time_ns')

def find_nearest(v1, v2):
    i1 = []
    i2 = []
    i = 0
    for v in v1:
        idx = (numpy.abs(v2-v).argmin())
        if idx not in i2:
            i2.append(idx)
            i1.append(i)
        i += 1
    return i1, i2
        

pylab.figure(1)
p = pylab.scatter(pend_time, psnr, label="PYCBC", color="blue", alpha=0.5, marker='x')
l = pylab.scatter(lend_time, lsnr, label="LAL", color="red", alpha=0.5, marker='+')
pylab.axhline(5.5)
for l in numpy.arange(64, 2048, 128):
    pylab.axvline(l)
pylab.title("SNR triggers over time")
pylab.xlabel("Time(s)")
pylab.ylabel("SNR")
pylab.legend()
pylab.savefig("snr_over_time.png")

pti, lti = find_nearest(pend_time, lend_time)
#pti = numpy.argsort(pend_time)
#lti = numpy.argsort(lend_time)
rel = (psnr[pti] - lsnr[lti]) / lsnr[lti]
relsig = (psig[pti]**0.5 - lsig[lti]**0.5) / lsig[lti]**0.5
cs = pmc[pti]


pylab.figure(4)
l = pylab.scatter( lend_time[lti], rel, label="LAL", c=cs, alpha=0.5)
pylab.colorbar()
for l in numpy.arange(64, 2048, 128):
    pylab.axvline(l)
pylab.title("relative difference  of SNR triggers over time")
pylab.xlabel("Time(s)")
pylab.ylabel("SNR")
pylab.savefig("snr_over_time_rel_mchirp.png")

pylab.figure(14)
l = pylab.scatter( lend_time[lti], rel, label="LAL", c=psnr[pti], alpha=0.5)
pylab.colorbar()
for l in numpy.arange(64, 2048, 128):
    pylab.axvline(l)
pylab.title("relative difference  of SNR triggers over time")
pylab.xlabel("Time(s)")
pylab.ylabel("SNR")
pylab.savefig("snr_over_time_rel_snr.png")

pylab.figure(5)
l = pylab.scatter(lend_time[lti], relsig, label="LAL", c=cs, alpha=0.5)
pylab.colorbar()
for l in numpy.arange(64, 2048, 128):
    pylab.axvline(l)
pylab.title("relative difference in sigma of triggers over time")
pylab.xlabel("Time(s)")
pylab.ylabel("SNR")
pylab.savefig("sigma_over_time_rel.png")

pylab.figure(3)
p = pylab.scatter(pend_time, psig, label="PYCBC", color="blue", alpha=0.5, marker='x')
l = pylab.scatter(lend_time, lsig, label="LAL", color="red", alpha=0.5, marker='+')
#l = pylab.scatter(lend_time, rel, label="LAL", c=cs, alpha=0.5)
#pylab.colorbar()
for l in numpy.arange(64, 2048, 128):
    pylab.axvline(l)
pylab.title("sigmasq of triggers over time")
pylab.xlabel("Time(s)")
pylab.ylabel("sigmasq")
pylab.legend()
pylab.savefig("sigma_over_time.png")

pylab.figure(2)
pylab.scatter(psnr, pchi, label="PYCBC", color="blue", alpha=0.5)
pylab.scatter(lsnr, lchi, label="LAL", color="red", alpha=0.5)
pylab.xlabel("snr")
pylab.ylabel("chisq")
pylab.legend()
pylab.savefig("snrchi.png")

pylab.show()

