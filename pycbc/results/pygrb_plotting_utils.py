import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import rc
import sys
import os
import argparse
import copy
import numpy
from pycbc.results import save_fig_with_metadata
#TODO: get rid of this (or move get_bestnr from pylal)
from pylal.coh_PTF_pyutils import get_bestnr, get_det_response

# Get rcParams
rc('font', size=14)

ptfcolormap = plt.cm.spring
ptfcolormap.set_over('g')

# =============================================================================
# Parse command line
# =============================================================================

def pygrb_plot_two_snrs_opts_parser(description='', version=None):
    parser = argparse.ArgumentParser(usage='', description=description)
    
    parser.add_argument("--version", action="version", version=version)
    
    parser.add_argument("-v", "--verbose", default=False, action="store_true",\
                      help="verbose output")
    
    parser.add_argument("-t", "--trig-file", action="store", \
                      default=None, required=True,\
                      help="The location of the trigger file")
    
    parser.add_argument("-I", "--inj-file", action="store", \
                      default=None, help="The location of the injection file")
    
    parser.add_argument("-a", "--segment-dir", action="store", 
                      required=True, help="directory holding buffer, on and "+\
                           "off source segment files.")
    
    parser.add_argument("-o", "--output-file", default=None, required=True,
                      help="Output file.")
    
    parser.add_argument("-O", "--zoomed-output-file", default=None, required=True,
                      help="Output file for a zoomed in verion of the plot.")
    
    parser.add_argument("-Q", "--chisq-index", action="store", type=float,\
                      default=4.0, help="chisq_index for newSNR calculation")
    
    parser.add_argument("-N", "--chisq-nhigh", action="store", type=float,\
                      default=3.0, help="nhigh for newSNR calculation")
    
    parser.add_argument("-B", "--sngl-snr-threshold", action="store", type=float,\
                      default=4.0, help="Single detector SNR threshold, the"+\
                      "two most sensitive detectors should have SNR above this")
    
    parser.add_argument("-d", "--snr-threshold", action="store", type=float,\
                      default=6.0, help="SNR threshold for recording triggers")
    
    parser.add_argument("-c", "--newsnr-threshold", action="store", type=float,\
                      default=None, help="NewSNR threshold for calculating the "+\
                      "chisq of triggers (based on value of auto and bank chisq"+\
                      " values. By default will take the same value as "+\
                      "snr-threshold")
    
    parser.add_argument("-A", "--null-snr-threshold", action="store",\
                      default="4.25,6",\
                      help="comma separated lower,higher null SNR thresholds, "+\
                           " for null SNR cut")
    
    parser.add_argument("-C", "--null-grad-thresh", action="store", type=float,\
                      default=20., help="Threshold above which to increase the,"+\
                      "values of the null SNR cut.")
    
    parser.add_argument("-D", "--null-grad-val", action="store", type=float,\
                      default=0.2, help="Rate the null SNR cut will increase"+\
                      "above the threshold.")
    
    parser.add_argument("-l", "--veto-directory",action="store",\
                      default=None,\
                      help="The location of the CATX veto files")
    
    parser.add_argument("-b", "--veto-category",action="store", type=int,\
                     default=None, help="Apply vetoes up to this level inclusive")
    
    return parser.parse_args()


# =============================================================================
# Format single detector chi-square data as numpy array and floor at 0.005
# =============================================================================

def format_single_chisqs(trigIfoCS, trigs, ifos):
  for ifo in ifos:
    trigIfoCS[ifo] = numpy.asarray(trigIfoCS[ifo])
    numpy.putmask(trigIfoCS[ifo], trigIfoCS[ifo]==0, 0.005)

  return trigIfoCS


# =============================================================================
# Extract trigger/injection data produced by PyGRB
# =============================================================================
class pygrb_filter_output(object):
  def __init__(self, trigs_or_injs, ifos, columns, output_type, chisq_index, chisq_nhigh, null_thresh, snrThresh, snglSnrThresh, newSnrThresh, nullGradThresh, nullGradVal, verbose=False):
    if verbose:
      sys.stdout.write("Extracting data related to the %s just loaded...\n" % output_type)
    # Initialize all content of self
    self.Time = None
    self.SNR = numpy.array(None)
    self.ReweightedSNR = None
    self.NullSNR = None
    self.Nullstat = None
    self.TraceSNR = None
    self.ChiSquare = numpy.array(None)
    self.BankVeto = None
    self.AutoVeto = None
    self.CoincSNR = None
    self.IfoSNR = dict((ifo, None) for ifo in ifos)
    self.FirstSNR = None
    self.SecondSNR = None
    self.ThirdSNR = None
    self.IfobankCS = dict((ifo, None) for ifo in ifos)
    self.IfoautoCS = dict((ifo, None) for ifo in ifos)
    self.IfostanCS = dict((ifo, None) for ifo in ifos)
    self.RelAmp1 = None
    self.Norm3 = None
    self.RelAmp2 = None
    self.Inclination = None
    # Exctract data and fill in content of self
    if trigs_or_injs is not None:
      # Work out if using sngl chisqs
      ifoAtt = { 'G1':'g', 'H1':'h1', 'H2':'h2', 'L1':'l', 'V1':'v', 'T1':'t' } 
      i = ifoAtt[ifos[0]]
  
      self.sngl_chisq      = 'chisq_%s' % i in columns
      self.sngl_bank_chisq = 'bank_chisq_%s' % i in columns
      self.sngl_cont_chisq = 'cont_chisq_%s' % i in columns
      
      # Set basic data
      self.Time          = numpy.asarray(trigs_or_injs.get_end())
      self.SNR           = numpy.asarray(trigs_or_injs.get_column('snr'))
      self.ReweightedSNR = [get_bestnr(t,q=chisq_index, n=chisq_nhigh,\
                                       null_thresh=null_thresh,snr_threshold=snrThresh,\
                                       sngl_snr_threshold = snglSnrThresh,\
                                       chisq_threshold = newSnrThresh,\
                                       null_grad_thresh = nullGradThresh,\
                                       null_grad_val = nullGradVal) for t in trigs_or_injs]
      self.ReweightedSNR = numpy.array(self.ReweightedSNR)
      self.NullSNR       = numpy.asarray(trigs_or_injs.get_null_snr())
      self.Nullstat      = numpy.asarray(trigs_or_injs.get_column('null_statistic'))
      self.TraceSNR      = numpy.asarray(trigs_or_injs.get_column('null_stat_degen'))
    
      # Get chisq data
      self.ChiSquare = numpy.asarray(trigs_or_injs.get_column('chisq'))
      self.BankVeto  = numpy.asarray(trigs_or_injs.get_column('bank_chisq')) 
      self.AutoVeto  = numpy.asarray(trigs_or_injs.get_column('cont_chisq'))
      numpy.putmask(self.ChiSquare, self.ChiSquare==0, 0.005)
      numpy.putmask(self.BankVeto, self.BankVeto==0, 0.005)
      numpy.putmask(self.AutoVeto, self.AutoVeto==0, 0.005)
       
      # Get single detector data
      self.CoincSNR  = (trigs_or_injs.get_column('coinc_snr'))
      self.IfoSNR    = dict((ifo, trigs_or_injs.get_sngl_snr(ifo)) for ifo in ifos)
      tmp           = numpy.sort(self.IfoSNR.values(), 0)
      if len(ifos) > 0:
        self.FirstSNR  = tmp[-1,:]
      if len(ifos) > 1:
        self.SecondSNR = tmp[-2,:]
      if len(ifos) > 2:
        self.ThirdSNR  = tmp[-3,:]
      if self.sngl_bank_chisq:
        self.IfobankCS = trigs_or_injs.get_sngl_bank_chisqs(ifos)
        self.IfobankCS = format_single_chisqs(self.IfobankCS, trigs_or_injs, ifos)
      if self.sngl_cont_chisq:
        self.IfoautoCS = trigs_or_injs.get_sngl_cont_chisqs(ifos)
        self.IfoautoCS = format_single_chisqs(self.IfoautoCS, trigs_or_injs, ifos)
      if self.sngl_chisq:
        self.IfostanCS = trigs_or_injs.get_sngl_chisqs(ifos)
        self.IfostanCS = format_single_chisqs(self.IfostanCS, trigs_or_injs, ifos)
      
      # Initiate amplitude generator
      numAmp  = 4
      amplitudes = xrange(1,numAmp+1)
    
      # Get amplitude terms
      Amp              = dict((amp,\
                               numpy.asarray(trigs_or_injs.get_column('amp_term_%d' % amp)))\
                               for amp in amplitudes)
      #print numpy.count_nonzero(Amp[1]) # All trigAmp's are 0
      #print numpy.count_nonzero(Amp[2]) # Hence the 3 warnings
      #print numpy.count_nonzero(Amp[3])
      #print numpy.count_nonzero(Amp[4])
      self.RelAmp1     = numpy.sqrt((Amp[1]**2 + Amp[2]**2)/\
                                    (Amp[3]**2 + Amp[4]**2))
      GammaR           = Amp[1] - Amp[4]
      GammaI           = Amp[2] + Amp[3]
      DeltaR           = Amp[1] + Amp[4]
      DeltaI           = Amp[3] - Amp[2]
      Norm1            = DeltaR*DeltaR + DeltaI*DeltaI
      Norm2            = GammaR*GammaR + GammaI*GammaI
      self.Norm3       = ((Norm1**0.25) + (Norm2**0.25))**2
      AmpPlus          = (Norm1)**0.5 + (Norm2)**0.5
      AmpCross         = abs((Norm1)**0.5 - (Norm2)**0.5)
      self.RelAmp2     = AmpPlus/AmpCross
      self.Inclination = AmpCross/self.Norm3
      
      num_trigs_or_injs = len(trigs_or_injs)
      if num_trigs_or_injs<1:
        sys.stderr.write("WARNING: No %s found." % output_type)
      elif num_trigs_or_injs >= 1 and verbose:
        sys.stdout.write("%d %s foundd.\n"\
                         % (num_trigs_or_injs, output_type))
      if(output_type == "triggers"):
        Sigma  = trigs_or_injs.get_sigmasqs()
        self.SigmaTot = numpy.zeros(num_trigs_or_injs)
        # Get antenna response based parameters
        self.Longitude = numpy.degrees(trigs_or_injs.get_column('ra'))
        self.Latitude  = numpy.degrees(trigs_or_injs.get_column('dec'))
        self.fResp = dict((ifo, numpy.empty(num_trigs_or_injs)) for ifo in ifos)
        for i in xrange(num_trigs_or_injs):
          # Calculate fResp for each IFO is we haven't done so already
          fPlus,fCross     = get_det_response(self.Longitude[i], self.Latitude[i],\
                                              self.Time[i])
          for ifo in ifos:
            self.fResp[ifo][i]   = sum(numpy.array([fPlus[ifo], fCross[ifo]])**2)
            self.SigmaTot[i]    += Sigma[ifo][i] * self.fResp[ifo][i]
        
        for ifo in ifos:
          self.fResp[ifo] = self.fResp[ifo].mean()
        
        # Normalise trigSigma
        self.SigmaTot = numpy.array(self.SigmaTot)
        for ifo in ifos:
          Sigma[ifo] = numpy.asarray(Sigma[ifo]) / self.SigmaTot
        
        self.SigmaMean = {}
        self.SigmaMax  = {}
        self.SigmaMin  = {}
        for ifo in ifos:
          try:
            self.SigmaMean[ifo] = Sigma[ifo].mean()
            self.SigmaMax[ifo]  = Sigma[ifo].max()
            self.SigmaMin[ifo]  = Sigma[ifo].min()
          except ValueError:
            self.SigmaMean[ifo] = 0
            self.SigmaMax[ifo]  = 0
            self.SigmaMin[ifo]  = 0
    if verbose:
      sys.stdout.write("%s parameters extractedd\n"\
                       % (output_type))

# =============================================================================
# Plot contours in a plot were SNR is on the horizontal axis
# =============================================================================

def contour_plotter(axis, snr_vals, contours, colors, vert_spike=False):
  for i in range(len(contours)):
    plot_vals_x = []
    plot_vals_y = []
    if vert_spike:
      for j in range(len(snr_vals)):
        if contours[i][j] > 1E-15:
          # Workaround to ensure that the vertical spike is shown on veto plots
          if not plot_vals_x:
            plot_vals_x.append(snr_vals[j])
            plot_vals_y.append(0.1)
          plot_vals_x.append(snr_vals[j])
          plot_vals_y.append(contours[i][j])
    else:
      plot_vals_x = snr_vals
      plot_vals_y = contours[i]
    axis.plot(plot_vals_x,plot_vals_y,colors[i])


# =============================================================================
# Master plotting funtion: fits all plotting needs in this script
# =============================================================================

def pygrb_plotter(trig_x, trig_y, inj_x, inj_y, injFile, xlabel, ylabel, \
                  fig_path, snr_vals=None, conts=None, \
                  shade_cont_value=None, colors=None, vert_spike=False,
                  xlims=None, ylims=None, use_logs=True, verbose=False):
  if verbose:
    fig_name = os.path.split(os.path.abspath(fig_path))[1]
    sys.stdout.write(" * %s (%s vs %s)...\n" % (fig_name, xlabel, ylabel))
  # Plot trigger-related quantities
  fig = plt.figure()
  ax = fig.gca()
  if use_logs:
    ax.loglog(trig_x, trig_y, 'bx')
  else:
    ax.plot(trig_x, trig_y, 'bx')
  ax.grid()
  # Plot injection-related quantities
  if injFile: 
    if use_logs:
      ax.loglog(inj_x, inj_y, 'r+')
    else:
      ax.plot(inj_x, inj_y, 'r+')
  # Plot contours
  if conts is not None:
    contour_plotter(ax, snr_vals, conts, colors, vert_spike=vert_spike)
  # Add shading above a specific contour (typically used for vetoed area)
  if shade_cont_value is not None:
    limy = ax.get_ylim()[1]
    polyx = copy.deepcopy(snr_vals)
    polyy = copy.deepcopy(conts[shade_cont_value])
    polyx = numpy.append(polyx,[max(snr_vals), min(snr_vals)])
    polyy = numpy.append(polyy, [limy, limy])
    ax.fill(polyx, polyy, color='#dddddd')
  # Axes: labels and limits 
  ax.set_xlabel(xlabel)
  ax.set_ylabel(ylabel)
  if xlims:
    ax.set_xlim(xlims)
  if ylims:
    ax.set_ylim(ylims)
  # Wrap up
  save_fig_with_metadata(fig, fig_path)
                         #fig_kwds=fig_kwds, title=plot_title,
                         #cmd=' '.join(sys.argv),
                         #caption=plot_caption)
  plt.close()
  return


# TODO: delete this when the pygrb_pp_workflow is complete
# =============================================================================
# Master plotting funtion: fits all plotting needs in this script
# =============================================================================

def sbv_plot_handler(trig_x, trig_y, inj_x, inj_y, injFile, xlabel, ylabel, \
                     outdir, fig_tag, snr_vals=None, conts=None, \
                     shade_cont_value=None, colors=None, vert_spike=False,
                     xlims=None, ylims=None, use_logs=True, verbose=False):
  if injFile: 
    file_label = "inj"
  else:
    file_label = "noinj"
  if verbose:
    sys.stdout.write(" * %s_%s.png (%s vs %s)...\n" % (fig_tag, file_label, xlabel, ylabel))
  # Plot trigger-related quantities
  fig = plt.figure()
  ax = fig.gca()
  if use_logs:
    ax.loglog(trig_x, trig_y, 'bx')
  else:
    ax.plot(trig_x, trig_y, 'bx')
  ax.grid()
  # Plot injection-related quantities
  if injFile: 
    if use_logs:
      ax.loglog(inj_x, inj_y, 'r+')
    else:
      ax.plot(inj_x, inj_y, 'r+')
  # Plot contours
  if conts is not None:
    contour_plotter(ax, snr_vals, conts, colors, vert_spike=vert_spike)
  # Add shading above a specific contour (typically used for vetoed area)
  if shade_cont_value is not None:
    limy = ax.get_ylim()[1]
    polyx = copy.deepcopy(snr_vals)
    polyy = copy.deepcopy(conts[shade_cont_value])
    polyx = numpy.append(polyx,[max(snr_vals), min(snr_vals)])
    polyy = numpy.append(polyy, [limy, limy])
    ax.fill(polyx, polyy, color='#dddddd')
  # Axes: labels and limits 
  ax.set_xlabel(xlabel)
  ax.set_ylabel(ylabel)
  if xlims:
    ax.set_xlim(xlims)
  if ylims:
    ax.set_ylim(ylims)
  # Wrap up
  fig.savefig("%s/%s_%s.png" % (outdir, fig_tag, file_label), bbox_inches='tight')
  plt.close()
  return
