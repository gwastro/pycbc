import matplotlib.pyplot as plt
import sys
import copy
import numpy

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
