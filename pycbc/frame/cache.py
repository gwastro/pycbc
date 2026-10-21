# Copyright (C) 2026 Alex Nitz
#
# This program is free software; you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the
# Free Software Foundation; either version 3 of the License, or (at your
# option) any later version.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
# Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
"""
Lists of data files (caches) and the file formats that can be read from them.
"""

import glob
import re
from collections import namedtuple
from pathlib import Path
from urllib.parse import urlparse

from pycbc.frame import gwf, gwosc_hdf

# Each supported format is a module providing the functions
# is_format, span, channel_info and read, which act on a single file.
_FORMATS = (gwf, gwosc_hdf)


def file_format(path):
    """Return the format module that can read `path`, or None."""
    for fmt in _FORMATS:
        if fmt.is_format(path):
            return fmt
    return None


# File naming convention of LIGO-T050017: OBS-DESCRIPTION-START-DURATION.ext
_T050017 = re.compile(
    r"^(?P<obs>[^-]+)-(?P<desc>[^-]+)-(?P<start>\d+(\.\d+)?)"
    r"-(?P<duration>\d+(\.\d+)?)\.[^.]+$"
)


class CacheEntry(
    namedtuple("CacheEntry", "observatory description start duration path")
):
    """One data file and the span of time it contains."""

    __slots__ = ()

    @property
    def end(self):
        return self.start + self.duration

    @property
    def url(self):
        return "file://localhost" + self.path

    @classmethod
    def from_file(cls, path):
        """
        Describe a file, using its name if it follows LIGO-T050017 and
        its content otherwise.
        """
        path = str(Path(path).resolve())
        match = _T050017.match(Path(path).name)
        if match:
            return cls(
                match["obs"],
                match["desc"],
                float(match["start"]),
                float(match["duration"]),
                path,
            )
        fmt = file_format(path)
        if fmt is None:
            raise TypeError(f"{path} is not a supported data file")  # noqa: TRY003
        start, duration = fmt.span(path)
        return cls("-", "-", start, duration, path)

    @classmethod
    def from_lal_line(cls, line):
        """Parse a line of a LAL-format cache file."""
        obs, desc, start, duration, url = line.split()
        return cls(obs, desc, float(start), float(duration), urlparse(url).path)


class FrameCache(list):
    """A list of `CacheEntry`, sorted by start time."""

    def sieve(self, regex=None, start=None, end=None):
        """
        Return the entries whose URL matches `regex` and which overlap
        the interval [`start`, `end`).
        """
        pattern = re.compile(regex) if regex else None
        return FrameCache(
            e
            for e in self
            if (pattern is None or pattern.search(e.url))
            and (start is None or e.end > start)
            and (end is None or e.start < end)
        )


def locations_to_cache(locations, latest=False):
    """
    Return a cumulative cache file build from the list of locations

    Parameters
    ----------
    locations : list
        A list of strings containing files, globs, or cache files used to
        build a combined cache.
    latest : Optional, {False, Boolean}
        Only return a cache with the most recent frame in the locations.
        If false, all results are returned.

    Returns
    -------
    cache : FrameCache
        A cumulative cache containing the files derived from the list of
        locations, sorted by start time.

    """
    entries = {}
    for source in locations:
        flist = glob.glob(source)  # noqa: PTH207
        if latest:

            def relaxed_getctime(fn):
                # when building a cache from a directory of temporary
                # low-latency frames, files might disappear between
                # the glob() and getctime() calls
                try:
                    return Path(fn).stat().st_ctime
                except OSError:
                    return 0

            if not flist:
                raise ValueError("no frame or cache files found in " + source)
            flist = [max(flist, key=relaxed_getctime)]

        for file_path in flist:
            if Path(file_path).suffix in [".lcf", ".cache"]:
                with Path(file_path).open() as cache_file:
                    new = [
                        CacheEntry.from_lal_line(line)
                        for line in cache_file
                        if line.strip() and not line.startswith("#")
                    ]
            elif file_format(file_path) is not None:
                new = [CacheEntry.from_file(file_path)]
            else:
                raise TypeError("Invalid location name")  # noqa: TRY003
            # Copies of one file in different directories hold the same data,
            # so only the first copy found is used
            for e in new:
                entries.setdefault(Path(e.path).name, e)
    return FrameCache(sorted(entries.values(), key=lambda e: (e.start, e.path)))
