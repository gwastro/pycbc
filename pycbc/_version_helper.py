# Based on generateGitID.sh by Reinhard Prix
#
# Copyright (C) 2009,2010, Adam Mercer <adam.mercer@ligo.org>
# Copyright (C) 2009,2010, Nickolas Fotopoulos <nvf@gravity.phys.uwm.edu>
# Copyright (C) 2008,2009, John T. Whelan <john.whelan@ligo.org>
# Copyright (C) 2008, Reinhard Prix <reinhard.ligo.org>
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 3 of the License, or (at
# your option) any later version.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
# General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.

__author__ = 'Adam Mercer <adam.mercer@ligo.org>'

import exceptions
import os
import time
import subprocess
import re


class GitInfo(object):
    def __init__(self):
        date = None
        hash = None
        branch = None
        tag = None
        author = None
        committer = None
        status = None
        builder = None
        build_date = None


class GitInvocationError(exceptions.LookupError):
    pass


def call(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
         on_error='ignore', returncode=False):
    """Run the given command (with shell=False) and return the output as a
    string.

    Strips the output of enclosing whitespace.

    If the return code is non-zero, throw GitInvocationError.
    """
    # start external command process
    p = subprocess.Popen(command, stdout=subprocess.PIPE,
                         stderr=subprocess.PIPE)

    # get outputs
    out, _ = p.communicate()

    # throw exception if process failed
    if p.returncode != 0 and on_error == 'raise':
        raise GitInvocationError('Failed to run "%s"' % " ".join(command))

    if returncode:
        return out.strip(), p.returncode
    else:
        return out.strip()


def get_build_name(git_path='git'):
    """Find the username of the current builder
    """
    name,retcode = call(('git', 'config', 'user.name'), returncode=True)
    if retcode:
        name = "Unknown User"
    email,retcode = call(('git', 'config', 'user.email'), returncode=True)
    if retcode:
        email = ""
    return "%s <%s>" % (name, email)


def get_build_date():
    """Returns the current datetime as the git build date
    """
    return time.strftime('%Y-%m-%d %H:%M:%S +0000', time.gmtime())


def get_last_commit(git_path='git'):
    """Returns the details of the last git commit

    Returns a tuple (hash, date, author name, author e-mail,
    committer name, committer e-mail).
    """
    hash_, udate, aname, amail, cname, cmail = (
        call((git_path, 'log', '-1',
              '--pretty=format:%H,%ct,%an,%ae,%cn,%ce')).split(","))
    date = time.strftime('%Y-%m-%d %H:%M:%S +0000', time.gmtime(float(udate)))
    author = '%s <%s>' % (aname, amail)
    committer = '%s <%s>' % (cname, cmail)
    return hash_, date, author, committer


def get_git_branch(git_path='git'):
    """Returns the name of the current git branch
    """
    branch_match = call((git_path, 'rev-parse', '--symbolic-full-name', 'HEAD'))
    if branch_match == "HEAD":
        return None
    else:
        return os.path.basename(branch_match)


def get_git_tag(hash_, git_path='git'):
    """Returns the name of the current git tag
    """
    tag, status = call((git_path, 'describe', '--exact-match',
                        '--tags', hash_), returncode=True)
    if status == 0:
        return tag
    else:
        return None

def get_git_status(git_path='git'):
    """Returns the state of the git working copy
    """
    status_output = subprocess.call((git_path, 'diff-files', '--quiet'))
    if status_output != 0:
        return 'UNCLEAN: Modified working tree'
    else:
        # check index for changes
        status_output = subprocess.call((git_path, 'diff-index', '--cached',
                                         '--quiet', 'HEAD'))
        if status_output != 0:
            return 'UNCLEAN: Modified index'
        else:
            return 'CLEAN: All modifications committed'


def generate_git_version_info():
    """Query the git repository information to generate a version module.
    """
    info = GitInfo()
    git_path = call(('which', 'git'))

    # get build info
    info.builder = get_build_name()
    info.build_date = get_build_date()

    # parse git ID
    info.hash, info.date, info.author, info.committer = (
        get_last_commit(git_path))

    # determine branch
    info.branch = get_git_branch(git_path)

    # determine tag
    info.tag = get_git_tag(info.hash, git_path)

    # determine version
    if info.tag:
        info.version = info.tag.strip('v')
        info.release = not re.search('[a-z]', info.version.lower())
    else:
        info.version = info.hash[:6]
        info.release = False

    # refresh index
    call((git_path, 'update-index', '-q', '--refresh'))

    # check working copy for changes
    info.status = get_git_status(git_path)

    return info
