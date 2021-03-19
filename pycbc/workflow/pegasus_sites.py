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

def add_local_site(sitecat, local_path, local_url): 
    local = Site("local", arch=Arch.X86_64, os_type=OS.LINUX)

    local_dir = Directory(Directory.SHARED_SCRATCH,
                          path=os.path.join(local_path, 'local-site-scratch'))
    local_file_serv = FileServer(urljoin(local_url, 'local-site-scratch'),
                                 Operation.ALL)
    local_dir.add_file_servers(local_file_serv)
    local.add_directories(local_dir)

    local_dir = Directory(Directory.LOCAL_STORAGE, path=local_path)
    local_file_serv = FileServer(local_url, Operation.ALL)
    local.add_directories(local_dir)

    local.add_profiles(Namespace.PEGASUS, key="style", value="condor")
    local.add_profiles(Namespace.CONDOR, key="getenv", value="True")
    local.add_profiles(Namespace.CONDOR, key="should_transfer_files",
                       value="Yes")
    local.add_profiles(Namespace.CONDOR, key="when_to_transfer_output",
                       value="ON_EXIT_OR_EVICT")

def add_condorpool_site(sitecat, local_path, local_url):
    site = Site("condorpool", arch=Arch.X86_64, os_type=OS.LINUX)

    local_dir = Directory(Directory.SHARED_SCRATCH,
                          path=os.path.join(local_path, 'local-site-scratch'))
    local_file_serv = FileServer(urljoin(local_url, 'local-site-scratch'),
                                 Operation.ALL)
    local_dir.add_file_servers(local_file_serv)
    site.add_directories(local_dir)
    
    local_dir = Directory(Directory.LOCAL_STORAGE, path=local_path)
    local_file_serv = FileServer(local_url, Operation.ALL)
    site.add_directories(local_dir)

    site.add_profiles(Namespace.PEGASUS, key="style", value="condor")
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

def add_nonfsio_site(sitecat):
    site = Site("nonfsio", arch=Arch.X86_64, os_type=OS.LINUX)
    site.add_profiles(Namespace.PEGASUS, key="style", value="condor")
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

def add_osg_site(sitecat):
    site = Site("osg", arch=Arch.X86_64, os_type=OS.LINUX)
    site.add_profiles(Namespace.PEGASUS, key="style", value="condor")
    site.add_profiles(Namespace.CONDOR, key="should_transfer_files",
                      value="Yes")
    site.add_profiles(Namespace.CONDOR, key="when_to_transfer_output",
                      value="ON_EXIT_OR_EVICT")
    site.add_profiles(Namespace.CONDOR, key="+OpenScienceGrid",
                      value="True")

def add_site(sitecat, sitename):
    curr_dir = os.getcwd()
    local_url = urljoin('file:', pathname2url(curr_dir))
    if sitename == 'local':
        add_local_site(sitecat, curr_dir, local_url)   
    elif sitename == 'condorpool':
        add_condorpool_site(sitecat, curr_dir, local_url)
    elif sitename == 'nonfsio':
        add_nonfsio_site(sitecat)
    elif sitename == 'osg':
        add_osg_site(sitecat)
    else:
        raise ValueError("Do not recognize site {}".format(sitename))

