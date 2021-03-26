# Copyright (C) 2021 The PyCBC development team

#
# =============================================================================
#
#                                   Preamble
#
# =============================================================================
#
""" This module provides default site catalogs, which should be suitable for
most use cases. You can override individual details here. It should also be
possible to implement a new site, but not sure how that would work in practice.
"""

import os
from urllib.parse import urljoin
from urllib.request import pathname2url
from Pegasus.api import Directory, FileServer, Site, Operation, Namespace
from Pegasus.api import Arch, OS

def add_site_pegasus_profile(site, cp):
    # Add global profile information
    if cp.has_section('pegasus_profile'):
        add_ini_site_profile(site, cp, 'pegasus_profile')
    # Add site-specific profile information
    if cp.has_section('pegasus_profile-{}'.format(site.name)):
        add_ini_site_profile(site, cp, 'pegasus_profile-{}'.format(site.name))

def add_ini_site_profile(site, cp, sec):
    for opt in cp.options(sec):
        namespace = opt.split('|')[0]
        if namespace == 'pycbc' or namespace == 'container':
            continue

        value = cp.get(sec, opt).strip()
        key = opt.split('|')[1]
        site.add_profiles(Namespace(namespace), key=key, value=value)

def add_local_site(sitecat, cp, local_path, local_url): 
    # local_url must end with a '/'
    if not local_url.endswith('/'):
        local_url = local_url + '/'

    local = Site("local", arch=Arch.X86_64, os_type=OS.LINUX)
    add_site_pegasus_profile(local, cp)

    local_dir = Directory(Directory.SHARED_SCRATCH,
                          path=os.path.join(local_path, 'local-site-scratch'))
    local_file_serv = FileServer(urljoin(local_url, 'local-site-scratch'),
                                 Operation.ALL)
    local_dir.add_file_servers(local_file_serv)
    local.add_directories(local_dir)

    local_dir = Directory(Directory.LOCAL_STORAGE, path=local_path)
    local_file_serv = FileServer(local_url, Operation.ALL)
    local_dir.add_file_servers(local_file_serv)
    local.add_directories(local_dir)

    local.add_profiles(Namespace.PEGASUS, key="style", value="condor")
    local.add_profiles(Namespace.PEGASUS, key="transfer.force",
                       value="true")
    local.add_profiles(Namespace.PEGASUS, key="transfer.links",
                       value="true")
    local.add_profiles(Namespace.CONDOR, key="getenv", value="True")
    local.add_profiles(Namespace.CONDOR, key="should_transfer_files",
                       value="Yes")
    local.add_profiles(Namespace.CONDOR, key="when_to_transfer_output",
                       value="ON_EXIT_OR_EVICT")
    sitecat.add_sites(local)

def add_condorpool_site(sitecat, cp, local_path, local_url):
    # local_url must end with a '/'
    if not local_url.endswith('/'):
        local_url = local_url + '/'

    site = Site("condorpool", arch=Arch.X86_64, os_type=OS.LINUX)
    add_site_pegasus_profile(site, cp)

    #local_dir = Directory(Directory.SHARED_SCRATCH,
    #                      path=os.path.join(local_path,
    #                                        'condorpool-site-scratch'))
    #local_file_serv = FileServer(urljoin(local_url, 'condorpool-site-scratch'),
    #                             Operation.ALL)
    #local_dir.add_file_servers(local_file_serv)
    #site.add_directories(local_dir)
    
    #local_dir = Directory(Directory.LOCAL_STORAGE, path=local_path)
    #local_file_serv = FileServer(local_url, Operation.ALL)
    #local_dir.add_file_servers(local_file_serv)
    #site.add_directories(local_dir)

    site.add_profiles(Namespace.PEGASUS, key="style", value="condor")
    site.add_profiles(Namespace.PEGASUS, key="transfer.links",
                      value="true")
    site.add_profiles(Namespace.PEGASUS, key="data.configuration",
                      value="nonsharedfs")
    site.add_profiles(Namespace.PEGASUS, key="auxillary.local",
                      value="True")
    site.add_profiles(Namepsace.PEGASUS, key='transfer.bypass.input.staging',
                      value="True")
    site.add_profiles(Namespace.CONDOR, key="should_transfer_files", 
                      value="Yes")
    site.add_profiles(Namespace.CONDOR, key="when_to_transfer_output", 
                      value="ON_EXIT_OR_EVICT")
    site.add_profiles(Namespace.CONDOR, key="getenv", value="True")
    site.add_profiles(Namespace.CONDOR, key="+DESIRED_Sites", 
                      value="'nogrid'")
    site.add_profiles(Namespace.CONDOR, key="+IS_GLIDEIN", 
                      value="'False'")
    site.add_profiles(Namespace.CONDOR, key="+flock_local", 
                      value="True")
    sitecat.add_sites(site)

def add_nonfsio_site(sitecat, cp):
    site = Site("nonfsio", arch=Arch.X86_64, os_type=OS.LINUX)
    add_site_pegasus_profile(site, cp)
    site.add_profiles(Namespace.PEGASUS, key="style", value="condor")
    # FIXME: condorio or nonsharedfs?
    site.add_profiles(Namespace.PEGASUS, key="data.configuration",
                      value="condorio")
    site.add_profiles(Namespace.CONDOR, key="should_transfer_files",
                      value="Yes")
    site.add_profiles(Namespace.CONDOR, key="when_to_transfer_output",
                      value="ON_EXIT_OR_EVICT")
    site.add_profiles(Namespace.CONDOR, key="+DESIRED_Sites",
                      value="'nogrid'")
    site.add_profiles(Namespace.CONDOR, key="+IS_GLIDEIN",
                      value="'False'")
    site.add_profiles(Namespace.CONDOR, key="+flock_local",
                      value="True")
    sitecat.add_sites(site)

def add_osg_site(sitecat, cp):
    site = Site("osg", arch=Arch.X86_64, os_type=OS.LINUX)
    add_site_pegasus_profile(site, cp)
    site.add_profiles(Namespace.PEGASUS, key="style", value="condor")
    # FIXME: condorio or nonsharedfs? Here frame files must *not* be copied
    site.add_profiles(Namespace.PEGASUS, key="data.configuration",
                      value="nonsharedfs")
    site.add_profiles(Namespace.CONDOR, key="should_transfer_files",
                      value="Yes")
    site.add_profiles(Namespace.CONDOR, key="when_to_transfer_output",
                      value="ON_EXIT_OR_EVICT")
    site.add_profiles(Namespace.CONDOR, key="+OpenScienceGrid",
                      value="True")
    # On OSG failure rate is high
    site.add_profiles(Namespace.DAGMAN, key="retry", value="4")
    sitecat.add_sites(site)

def add_site(sitecat, sitename, cp, out_dir=None):
    if out_dir is None:
        out_dir = os.getcwd()
    local_url = urljoin('file:', pathname2url(out_dir))
    if sitename == 'local':
        add_local_site(sitecat, cp, out_dir, local_url)   
    elif sitename == 'condorpool':
        add_condorpool_site(sitecat, cp, out_dir, local_url)
    elif sitename == 'nonfsio':
        add_nonfsio_site(sitecat, cp)
    elif sitename == 'osg':
        add_osg_site(sitecat, cp)
    else:
        raise ValueError("Do not recognize site {}".format(sitename))

