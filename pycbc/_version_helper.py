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

import os
import time
import subprocess
import re
import distutils.version


class GitInfo(object):
    def __init__(self):
        self.date = None
        self.hash = None
        self.branch = None
        self.tag = None
        self.author = None
        self.committer = None
        self.status = None
        self.builder = None
        self.build_date = None


class GitInvocationError(LookupError):
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

    out = out.decode('utf-8')

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

def get_num_commits():
    return call(('git', 'rev-list', '--count', 'HEAD'))

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

def determine_latest_release_version():
    """Query the git repository for the last released version of the code.
    """
    git_path = call(('which', 'git'))

    # Get all tags
    tag_list = call((git_path, 'tag')).split('\n')

    # Reduce to only versions
    tag_list = [t[1:] for t in tag_list if t.startswith('v')]

    # Determine if indeed a tag and store largest
    latest_version = None
    latest_version_string = None
    re_magic = re.compile("\d+\.\d+\.\d+$")
    for tag in tag_list:
        # Is this a version string
        if re_magic.match(tag):
            curr_version = distutils.version.StrictVersion(tag)
            if latest_version is None or curr_version > latest_version:
                latest_version = curr_version
                latest_version_string = tag

    return latest_version_string

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
        info.version = '0.0a' + get_num_commits()
        info.release = False

    # Determine *last* stable release
    info.last_release = determine_latest_release_version()

    # refresh index
    call((git_path, 'update-index', '-q', '--refresh'))

    # check working copy for changes
    info.status = get_git_status(git_path)

    return info
