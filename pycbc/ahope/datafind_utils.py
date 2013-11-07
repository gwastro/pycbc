import os,sys,optparse
from glue import datafind
from glue import segments,segmentsUtils,git_version


def run_datafind_local(cp):
    """
    Add documentation.
    """
    # First job is to do setup for the datafind jobs
    # First get the server name
    if 0:
        # Placeholder for getting server name from the config file
        pass
    else:
        # Get the server name from the environment
        if os.environ.has_key("LIGO_DATAFIND_SERVER"):
            datafindServer = os.environ["LIGO_DATAFIND_SERVER"]
        else:
            errMsg = "Trying to obtain the ligo datafind server url from "
            errMsg += "the environment, ${LIGO_DATAFIND_SERVER}, but that "
            errMsg += "variable is not populated."
            raise ValueError(errMsg)

    # verify authentication options
    if not datafindServer.endswith("80"):
        cert_file, key_file = datafind.find_credential()
    else:
        cert_file, key_file = None, None

    # Is a port specified in the server URL
    server, port = datafindServer.split(':',1)
    if port == "":
        port = None
    else:
        port = int(port)

    # Open connection to the datafind server
    if cert_file and key_file:
        #HTTPS connection
        connection =\
            datafind.GWDataFindHTTPSConnection(host=server, port=port, \
                                   cert_file=cert_file, key_file=key_file)
    else:
        # HTTP connection
        connection =\
            datafind.GWDataFindHTTPConnection(host=server, port=port)

    # Now ready to loop over the input segments
    ifos = scienceSegmentsAllIfos.keys()
    for ifo, scienceSegsIfo in scienceSegmentsAllIfos.items():
        observatory = ifo[0].upper()
        # Get the type from the config file
        type = 'BLOOPY'
        for seg in scienceSegsIfo:
            # FIXME: Do these need to be integers?
            startTime = seg[0]
            endTime = seg[1]
        # FIXME: Take KWargs from config
        try:
            dfCache = connection.find_frame_urls(observatory, type, start_time,\
                             end_time, urltype="file", on_gaps="error")
        except RuntimeError:
            raise
        # cache.tofile(sys.stderr)


    
