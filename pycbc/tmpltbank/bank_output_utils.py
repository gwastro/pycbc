from __future__ import division
import numpy
from lal import LAL_PI, LAL_MTSUN_SI, LAL_TWOPI, LAL_GAMMA
from glue.ligolw import ligolw
from glue.ligolw import table
from glue.ligolw import lsctables
from glue.ligolw import ilwd
from glue.ligolw import utils as ligolw_utils
from glue.ligolw.utils import process as ligolw_process
from pycbc import pnutils
from pycbc.tmpltbank.lambda_mapping import *

def return_empty_sngl():
    """
    Function to create a SnglInspiral object where all columns are populated
    but all are set to values that test False (ie. strings to '', floats/ints
    to 0, ...). This avoids errors when you try to create a table containing
    columns you don't care about, but which still need populating. NOTE: This
    will also produce a process_id and event_id with 0 values. For most
    applications these should be set to their correct values.

    Returns
    --------
    lsctables.SnglInspiral
        The "empty" SnglInspiral object.
    """

    sngl = lsctables.SnglInspiral()
    cols = lsctables.SnglInspiralTable.validcolumns
    for entry in cols.keys():
        if cols[entry] in ['real_4','real_8']:
            setattr(sngl,entry,0.)
        elif cols[entry] == 'int_4s':
            setattr(sngl,entry,0)
        elif cols[entry] == 'lstring':
            setattr(sngl,entry,'')
        elif entry == 'process_id':
            sngl.process_id = ilwd.ilwdchar("sngl_inspiral:process_id:0")
        elif entry == 'event_id':
            sngl.event_id = ilwd.ilwdchar("sngl_inspiral:event_id:0")
        else:
            raise ValueError("Column %s not recognized" %(entry) )
    return sngl

def convert_to_sngl_inspiral_table(params, proc_id):
    '''
    Convert a list of m1,m2,spin1z,spin2z into a sngl_inspiral format template
    bank.

    Parameters
    -----------
    params : list of tuples
        Each entry in the params list should be a tuple/list/numpy array
        whose entries contain [mass1,mass2,spin1z,spin2z]
    proc_id : ilwd char
        Process ID to add to each row of the sngl_inspiral table

    Returns
    ----------
    SnglInspiralTable
        The bank of templates in SnglInspiralTable format.
    '''
    sngl_inspiral_table = lsctables.New(lsctables.SnglInspiralTable)
    col_names = ['mass1','mass2','spin1z','spin2z']

    for values in params:
        tmplt = return_empty_sngl()

        tmplt.process_id = proc_id
        index = 0
        for value in values[0:4]:
            setattr(tmplt,col_names[index],value)
            index += 1
        tmplt.mtotal = tmplt.mass1 + tmplt.mass2
        tmplt.eta = tmplt.mass1 * tmplt.mass2 / (tmplt.mtotal * tmplt.mtotal)
        tmplt.mchirp = tmplt.mtotal * tmplt.eta**(3./5.)
        # Currently using ISCO frequency for termination
        tmplt.f_final = pnutils.f_SchwarzISCO(tmplt.mtotal)
        tmplt.template_duration = 0 # FIXME
        tmplt.event_id = sngl_inspiral_table.get_next_id()
        # FIXME: Add gamma values
        sngl_inspiral_table.append(tmplt)

    return sngl_inspiral_table

def calculate_ethinca_metric_comps(metricParams, ethincaParams, mass1, mass2, 
                                   spin1z=0., spin2z=0.):
    """
    Calculate the Gamma components needed to use the ethinca metric.
    At present this outputs the standard TaylorF2 metric over the end time 
    and chirp times \tau_0 and \tau_3.
    A desirable upgrade might be to use the \chi coordinates [defined WHERE?] 
    for metric distance instead of \tau_0 and \tau_3.
    The lower frequency cutoff is currently hard-coded to be the same as the 
    bank layout options fLow and f0 (which must be the same as each other).

    Parameters
    -----------
    metricParams : metricParameters instance
        Structure holding all the options for construction of the metric
        and the eigenvalues, eigenvectors and covariance matrix
        needed to manipulate the space.
    ethincaParams: ethincaParameters instance
        Structure holding options relevant to the ethinca metric computation.

    Returns
    --------
    numpy_array
        A numpy array holding 6 independent metric components in 
        (end_time, tau_0, tau_3) coordinates to be stored in the Gamma0-5 
        slots of a SnglInspiral object.
    """
    if (float(spin1z) != 0. or float(spin2z) != 0.):
        raise NotImplementedError("Ethinca cannot at present be calculated "
                                  "for nonzero component spins!")
    f0 = metricParams.f0
    if f0 != metricParams.fLow: 
        raise ValueError("If calculating ethinca the bank f0 value must be "
                         "equal to f-low!")
    if ethincaParams.fLow is not None and (
        ethincaParams.fLow != metricParams.fLow):
        raise NotImplementedError("An ethinca metric f-low different from the"
                                  " bank metric f-low is not supported!")

    twicePNOrder = ethinca_order_from_string(ethincaParams.pnOrder)

    piFl = LAL_PI * f0
    totalMass, eta = pnutils.mass1_mass2_to_mtotal_eta(mass1, mass2)
    totalMass = totalMass * LAL_MTSUN_SI
    v0cube = totalMass*piFl
    v0 = v0cube**(1./3.)

    # Get theoretical cutoff frequency and work out the closest 
    # frequency for which moments were calculated
    fMax_theor = pnutils.frequency_cutoff_from_name(
        ethincaParams.cutoff, mass1, mass2, spin1z, spin2z)
    fMaxes = metricParams.moments['J4'].keys()
    fMaxIdx = abs(numpy.array(fMaxes,dtype=float) - fMax_theor).argmin()
    fMax = fMaxes[fMaxIdx]

    # 3pN is a mess, so split it into pieces
    a0 = 11583231236531/200286535680 - 5*LAL_PI*LAL_PI - 107*LAL_GAMMA/14
    a1 = (-15737765635/130056192 + 2255*LAL_PI*LAL_PI/512)*eta
    a2 = (76055/73728)*eta*eta
    a3 = (-127825/55296)*eta*eta*eta
    alog = numpy.log(4*v0) # Log terms are tricky - be careful

    # Get the Psi coefficients
    Psi = [{},{}] #Psi = numpy.zeros([2,8,2],dtype=float)
    Psi[0][0,0] = 3/5
    Psi[0][2,0] = (743/756 + 11*eta/3)*v0*v0
    Psi[0][3,0] = 0.
    Psi[0][4,0] = (-3058673/508032 + 5429*eta/504 + 617*eta*eta/24)\
                    *v0cube*v0
    Psi[0][5,1] = (-7729*LAL_PI/126)*v0cube*v0*v0/3
    Psi[0][6,0] = (128/15)*(-3*a0 - a1 + a2 + 3*a3 + 107*(1+3*alog)/14)\
                    *v0cube*v0cube
    Psi[0][6,1] = (6848/35)*v0cube*v0cube/3
    Psi[0][7,0] = (-15419335/63504 - 75703*eta/756)*LAL_PI*v0cube*v0cube*v0

    Psi[1][0,0] = 0.
    Psi[1][2,0] = (3715/12096 - 55*eta/96)/LAL_PI/v0;
    Psi[1][3,0] = -3/2
    Psi[1][4,0] = (15293365/4064256 - 27145*eta/16128 - 3085*eta*eta/384)\
                    *v0/LAL_PI
    Psi[1][5,1] = (193225/8064)*v0*v0/3
    Psi[1][6,0] = (4/LAL_PI)*(2*a0 + a1/3 - 4*a2/3 - 3*a3 -107*(1+6*alog)/42)\
                    *v0cube
    Psi[1][6,1] = (-428/LAL_PI/7)*v0cube/3
    Psi[1][7,0] = (77096675/1161216 + 378515*eta/24192 + 74045*eta*eta/8064)\
                    *v0cube*v0

    # Set the appropriate moments
    Js = numpy.zeros([18,3],dtype=float)
    for i in range(18):
        Js[i,0] = metricParams.moments['J%d'%(i)][fMax]
        Js[i,1] = metricParams.moments['log%d'%(i)][fMax]
        Js[i,2] = metricParams.moments['loglog%d'%(i)][fMax]

    # Calculate the g matrix
    PNterms = [(0,0),(2,0),(3,0),(4,0),(5,1),(6,0),(6,1),(7,0)]
    PNterms = [term for term in PNterms if term[0] <= twicePNOrder]
    two_pi_flower_sq = LAL_TWOPI * f0 * LAL_TWOPI * f0

    # Now can compute the gamma values
    gammaVals = numpy.zeros([6],dtype=float)
    gammaVals[0] = 0.5 * two_pi_flower_sq * \
                    ( Js[(1,0)] - (Js[(4,0)]*Js[(4,0)]) )

    for m in [0, 1]:
        for k in PNterms:
            gammaVals[1+m] += 0.5 * two_pi_flower_sq * Psi[m][k] * \
                                ( Js[(9-k[0],k[1])]
                                - Js[(12-k[0],k[1])] * Js[(4,0)] )

    g = numpy.zeros([2,2],dtype=float)
    for (m,n) in [(0,0),(0,1),(1,1)]:
        for k in PNterms:
            for l in PNterms:
                g[m,n] += Psi[m][k] * Psi[n][l] * \
                        ( Js[(17-k[0]-l[0], k[1]+l[1])]
                        - Js[(12-k[0],k[1])] * Js[(12-l[0],l[1])] )
        g[m,n] = 0.5 * two_pi_flower_sq * g[m,n]
        g[n,m] = g[m,n]

    gammaVals[3] = g[0,0]
    gammaVals[4] = g[0,1]
    gammaVals[5] = g[1,1]

    return gammaVals

def output_sngl_inspiral_table(outputFile, tempBank, metricParams,
                               ethincaParams, programName="", optDict = {}, 
                               outdoc=None, **kwargs):
    """
    Function that converts the information produced by the various pyCBC bank
    generation codes into a valid LIGOLW xml file, containing a sngl_inspiral
    table and output this to file.
 
    Parameters
    -----------
    outputFile : string
        Name of the file that the bank will be written to
    tempBank : list of tuples
        Each entry in the tempBank list should be a tuple/list/numpy array
        whose entries contain [mass1,mass2,spin1z,spin2z]
    metricParams : metricParameters instance
        Structure holding all the options for construction of the metric
        and the eigenvalues, eigenvectors and covariance matrix
        needed to manipulate the space.
    ethincaParams: ethincaParameters instance
        Structure holding options relevant to the ethinca metric computation.
        NOTE: The computation is currently only valid for non-spinning systems
        and uses the TaylorF2 approximant.
    programName (key-word-argument) : string
        Name of the executable that has been run
    optDict (key-word argument) : dictionary
        Dictionary of the command line arguments passed to the program
    outdoc (key-word argument) : ligolw xml document
        If given add template bank to this representation of a xml document and
        write to disk. If not given create a new document.
    kwargs : key-word arguments
        All other key word arguments will be passed directly to 
        ligolw_process.register_to_xmldoc
    """
    if outdoc is None:
        outdoc = ligolw.Document()
        outdoc.appendChild(ligolw.LIGO_LW())
    proc_id = ligolw_process.register_to_xmldoc(outdoc, programName, optDict,
                                                **kwargs).process_id
    sngl_inspiral_table = \
            convert_to_sngl_inspiral_table(tempBank, proc_id)
    # Calculate Gamma components if needed
    if ethincaParams is not None and (ethincaParams.doEthinca == True):
        for sngl in sngl_inspiral_table:
            # Set tau_0 and tau_3 values needed for the calculation of 
            # ethinca metric distances
            (sngl.tau0,sngl.tau3) = pnutils.mass1_mass2_to_tau0_tau3(
                  sngl.mass1, sngl.mass2, metricParams.f0) 
            GammaVals = calculate_ethinca_metric_comps(
                metricParams, ethincaParams, sngl.mass1, sngl.mass2,
                sngl.spin1z, sngl.spin2z)
            # assign the output Gamma0-5 values to the sngl object
            for i in range(6): setattr(sngl, "Gamma"+str(i), GammaVals[i])

    outdoc.childNodes[0].appendChild(sngl_inspiral_table)

    # write the xml doc to disk
    proctable = table.get_table(outdoc, lsctables.ProcessTable.tableName)
    ligolw_utils.write_filename(outdoc, outputFile, \
                                gz=outputFile.endswith('.gz'))
