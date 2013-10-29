from __future__ import division
import numpy
from lal import LAL_PI, LAL_MTSUN_SI, LAL_TWOPI, LAL_GAMMA
from glue.ligolw import ligolw
from glue.ligolw import table
from glue.ligolw import lsctables
from glue.ligolw import ilwd
from glue.ligolw import utils as ligolw_utils
from glue.ligolw.utils import process as ligolw_process

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
        tmplt.f_final = (1/6.)**(3./2.) / (LAL_PI * tmplt.mtotal * LAL_MTSUN_SI)
        tmplt.template_duration = 0 # FIXME
        # FIXME: Add gamma values
        sngl_inspiral_table.append(tmplt)

    return sngl_inspiral_table

def calculate_ethinca_metric_comps(temp_params, moments, f0):
    """
    Calculate the Gamma components that are needed to do an ethinca metric.
    Note that while the rest of the bank can use 3.5PN or R2F4, this will only
    use 2PN F2 metric (feel free to update this!) The reason for this is that
    it might be better to use the \chi coordinates for metric distance instead
    of \tau_0 and \tau_3.
    NOTE: Large parts of this are stolen shamelessly from lalinspiral

    Parameters
    -----------
    temp_params : Tuple
        A tuple containing (mass1,mass2). 
    moments : dictionary
        This dictionary contains all the moments at various cutoff frequencies.
        See the get_moments function for information on this.
    f0 : float
        The fiducial frequency used to scale the metric calculations. The value
        is irrelevant, but most stay the same.

    Returns
    --------
    numpy_array
        A numpy array holding the 6 values to go into the various Gamma entries
        in the sngl_inspiral table.
    """
    # Get \tau_0 - \tau_3 coordinates of template
    piFl = LAL_PI * f0
    totalMass = (temp_params[0] + temp_params[1])
    eta = temp_params[0]*temp_params[1] / (totalMass*totalMass)
    totalMass = totalMass * LAL_MTSUN_SI
    tau0 = 5.0/(256.0*eta*(totalMass**(5./3.))*(piFl**(8./3.)))
    tau3 = LAL_PI/(8.0*eta*(totalMass**(2./3.))*(piFl**(5./3.)))
    t1 = LAL_TWOPI * f0 * tau0
    t2 = LAL_TWOPI * f0 * tau3

    # Figure out appropriate f_max
    tempFMax = (1/6.)**(3./2.) / (LAL_PI * totalMass)
    fMaxes = moments['J4'].keys()
    fMaxIdx = abs(numpy.array(fMaxes,dtype=float) -  tempFMax).argmin()
    fMax = fMaxes[fMaxIdx]

    # Calculate the t0/t3 moment components
    t0t3Moms = {}
    t0t3Moms['a01'] = 3./5.
    t0t3Moms['a21'] = 11. * LAL_PI/12.
    t0t3Moms['a22'] = 743./2016. * (25./(2.*LAL_PI*LAL_PI))**(1./3.)
    t0t3Moms['a31'] = -3./2.
    t0t3Moms['a41'] = 617. * LAL_PI * LAL_PI / 384.
    t0t3Moms['a42'] = 5429./5376. * (25.*LAL_PI/2.)**(1./3.)
    t0t3Moms['a43'] = 1.5293365/1.0838016 * \
                      (5./(4.*LAL_PI*LAL_PI*LAL_PI*LAL_PI))**(1./3.)

    # Get the Psi coefficients
    Psi = numpy.zeros([2,5],dtype=float)
    Psi[0][0] = t0t3Moms['a01']
    Psi[0][1] = 0.
    Psi[0][2] = t0t3Moms['a21']/t2 + t0t3Moms['a22']/3. * \
                (t2 * t2 / (t1 * t1))**(1./3.)
    Psi[0][3] = 0.
    Psi[0][4] = t0t3Moms['a41']/(t2*t2) \
              + t0t3Moms['a42']/(3.* (t1*t1*t2)**(1./3.)) \
              - t0t3Moms['a43']/3. * t2 / t1 * (t2 / t1)**(1./3.)
    Psi[1][0] = 0.
    Psi[1][1] = 0.
    Psi[1][2] = - t0t3Moms['a21']*t1/(t2*t2) + \
                2. * t0t3Moms['a22']/3. * (t1/t2)**(1./3.)
    Psi[1][3] = t0t3Moms['a31']
    Psi[1][4] = - 2. * t0t3Moms['a41']*t1 / (t2*t2*t2) - \
                t0t3Moms['a42']/3. * (t1/(t2*t2*t2*t2))**(1./3.) + \
                4. * t0t3Moms['a43']/3. * (t2/t1)**(1./3.)

    # Set the appropriate moments
    Js = []
    for i in range(18):
        Js.append(moments['J%d'%(i)][fMax])
    Js.append(moments['J%d'%(-1)][fMax])
    # The log moments are in here as well and will be needed for higher order
    # terms.
  
    # Calculate the g matrix
    PNorder = 5
    g = numpy.zeros([2,2],dtype=float)
    two_pi_flower_sq = LAL_TWOPI * f0 * LAL_TWOPI * f0
    for (m,n) in [(0,0),(0,1),(1,0),(1,1)]:
        for k in range(PNorder):
            for l in range(PNorder):
                g[m,n] += Psi[m][k] * Psi[n][l] * (\
                        Js[17-k-l] - Js[12-k] * Js[12-l] \
                        - ( Js[9-k] - Js[4] * Js[12-k] ) \
                        * ( Js[9-l] - Js[4] * Js[12-l] ) \
                        / ( Js[1] - Js[4] * Js[4] ) )
        g[m,n] = 0.5 * two_pi_flower_sq * g[m,n]

    # Now can compute the gamma values
    gammaVals = numpy.zeros([6],dtype=float)
    gammaVals[0] = 0.5 * two_pi_flower_sq * \
                  ( Js[1] - (Js[4]*Js[4]) )
    gammaVals[1] = 0.5 * two_pi_flower_sq * \
                  ( Psi[0][0]*(Js[9] - (Js[4]*Js[12]) ) \
                  + Psi[0][2]*(Js[7] - (Js[4]*Js[10]) ) \
                  + Psi[0][4]*(Js[5] - (Js[4]*Js[8]) ))
    gammaVals[2] = 0.5 * two_pi_flower_sq * \
                  ( Psi[1][2]*(Js[7] - (Js[4]*Js[10]) ) \
                  + Psi[1][3]*(Js[6] - (Js[4]*Js[9]) ) \
                  + Psi[1][4]*(Js[5] - (Js[4]*Js[8]) ))
    gammaVals[3] = g[0,0] + gammaVals[1] * gammaVals[1] / gammaVals[0]
    gammaVals[4] = g[0,1] + gammaVals[1] * gammaVals[2] / gammaVals[0]
    gammaVals[5] = g[1,1] + gammaVals[2] * gammaVals[2] / gammaVals[0]

    return gammaVals

def output_sngl_inspiral_table(outputFile, tempBank, moments, f0,\
                               calculate_ethinca_comps=False,
                               programName="", optDict = {}, **kwargs):
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
    moments : Dictionary
        This must be the moments returned by the function 
        pycbc.tmpltbank.determine_eigen_directions
    f0 : float
        The value of f0 used in the moments calculation
    calculate_ethinca_comps (key-word-argument) : Boolean
        If set to True, this will calculate the ethinca metric components that
        are needed when doing lalapps/ligolw_thinca coincidence. NOTE: These
        moments are only valid for non-spinning systems and are currently only
        calculated to 2PN order.
    programName (key-word-argument) : string
        Name of the executable that has been run
    optDict (key-word argument) : dictionary
        Dictionary of the command line arguments that were passed to the program
    kwargs : key-word arguments
        All other key word arguments will be passed directly to 
        ligolw_process.register_to_xmldoc
    """
    outdoc = ligolw.Document()
    outdoc.appendChild(ligolw.LIGO_LW())
    proc_id = ligolw_process.register_to_xmldoc(outdoc, programName, optDict,\
                                                **kwargs).process_id
    sngl_inspiral_table = \
            convert_to_sngl_inspiral_table(tempBank, proc_id)
    # Calculate Gamma components if needed
    if calculate_ethinca_comps:
        for sngl in sngl_inspiral_table:
            temp_params = (sngl.mass1, sngl.mass2)
            GammaVals = calculate_ethinca_metric_comps(\
                        temp_params, moments, f0)
            sngl.Gamma0 = GammaVals[0]
            sngl.Gamma1 = GammaVals[1]
            sngl.Gamma2 = GammaVals[2]
            sngl.Gamma3 = GammaVals[3]
            sngl.Gamma4 = GammaVals[4]
            sngl.Gamma5 = GammaVals[5]

    outdoc.childNodes[0].appendChild(sngl_inspiral_table)

    # write the xml doc to disk
    proctable = table.get_table(outdoc, lsctables.ProcessTable.tableName)
    ligolw_utils.write_filename(outdoc, outputFile)
