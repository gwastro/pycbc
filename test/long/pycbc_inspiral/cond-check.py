import pylab
import numpy
from pycbc.frame import read_frame
from pycbc.types import load_timeseries

lal_out_file = "H1-INSPIRAL_lalsuite_FULL_DATA-968605000-2048.gwf"

# Take only the first 16 seconds to use as a small comparison
lal_raw = read_frame(lal_out_file, 'H1:LDAS-STRAIN_RAW', start_time=968605000, duration="16")
lal_resample = read_frame(lal_out_file, 'H1:LDAS-STRAIN_RAW_RESAMP', start_time=968605000, duration="16")
lal_conditioned = read_frame(lal_out_file, 'H1:LDAS-STRAIN_FILTER', start_time=968605000, duration="16")

# Take only the first 16 seconds to use as a small comparison
start_offset = 8 # This is the offest between the output that PyCBC gives and lalapps_inspiral 
                 # lalapps_inspiral does not output the 8 second padding in the frame file
pycbc_raw = load_timeseries("raw_pycbc.npy")[start_offset*16384:(start_offset+16)*16384].astype(lal_raw.dtype)
pycbc_resample = load_timeseries("resampled_pycbc.npy")[start_offset*4096:(start_offset+16)*4096].astype(lal_raw.dtype)

# the lalapps_inspiral gwf file misreportes the epoch for some of the streams
pycbc_conditioned = load_timeseries("conditioned_pycbc.npy")[start_offset*4096:(start_offset+16)*4096].astype(lal_raw.dtype)
pycbc_conditioned._epoch = lal_conditioned._epoch

print float(pycbc_conditioned.start_time)

def p(a, b, an, bn):
    t = numpy.arange(0, len(a), 1)
    mindif = 1e-7

    pylab.figure()
    pylab.scatter(t, a, label=an, color="red", marker='x')
    pylab.scatter(t, b, label=bn, color="blue", marker='+')
    pylab.legend()
    pylab.savefig("cond-plot" + an + '-' + bn + ".png")

    pylab.figure()
    pylab.scatter(t, a-b)
    ymin = min( min(a-b), -1 * mindif)
    ymax = max( max(a-b), 1 * mindif)
    pylab.ylim(ymin, ymax)
    pylab.title("absolute difference %s %s" % (an, bn))
    pylab.savefig("cond-abs_diff" + an + '-' + bn + ".png")

    pylab.figure()
    rdif = (a - b) / a
    pylab.scatter(t, rdif)
    ymin = min( min(rdif), -1 * mindif)
    ymax = max( max(rdif), 1 * mindif)
    pylab.ylim(ymin, ymax)
    pylab.title("fraction difference %s %s" % (an, bn))
    pylab.savefig("cond-frac_diff" + an + '-' + bn + ".png")
    
    pylab.figure()   
    bins=[-1e-5, -1e-6, -1e-7, -1e-8, 1e-8, 1e-8, 1e-7, 1e-6, 1e-5]   
    pylab.hist(rdif,  bins=bins)
    pylab.gca().set_xscale('symlog', linthreshx=1e-8)
    pylab.title("Fractional Difference %s - %s" % (an, bn))
    pylab.savefig("cond-hist-frac_diff" + an + '-' +  bn + ".png")
    
    pylab.figure()
    pylab.hist(a-b, bins=bins)
    pylab.gca().set_xscale('symlog', linthreshx=1e-8)
    pylab.title("Absolute Difference (%s - %s)" % (an, bn))
    pylab.savefig("cond-hist-abs_diff" + an + '-' + bn + ".png")

p(lal_raw, pycbc_raw, "lal-raw", "pycbc-raw")
p(lal_resample, pycbc_resample, "lal-resampled", "pycbc-resampled")
p(lal_conditioned, pycbc_conditioned, "lal-conditioned", "pycbc-conditioned")
