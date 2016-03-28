# Copyright (C) 2012  Alex Nitz
#
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


#
# =============================================================================
#
#                                   Preamble
#
# =============================================================================
#
"""This file contains backported functionality from future versions of libraries.
Mostly done with monkey-patching.
"""

def block_diag(*arrs):
    import numpy as np
    from numpy import atleast_1d, atleast_2d, array
    if arrs == ():
        arrs = ([],)
    arrs = [np.atleast_2d(a) for a in arrs]

    bad_args = [k for k in range(len(arrs)) if arrs[k].ndim > 2]
    if bad_args:
        raise ValueError("arguments in the following positions have dimension "
                            "greater than 2: %s" % bad_args)

    shapes = np.array([a.shape for a in arrs])
    out = np.zeros(np.sum(shapes, axis=0), dtype=arrs[0].dtype)

    r, c = 0, 0
    for i, (rr, cc) in enumerate(shapes):
        out[r:r + rr, c:c + cc] = arrs[i]
        r += rr
        c += cc
    return out

def zpk2sos(z, p, k, pairing='nearest'):
    """Stolen from scipy, please kill me and upgrade scipy...
    """

    import numpy as np
    from numpy import zeros
    from scipy.signal import zpk2tf, lfilter

    valid_pairings = ['nearest', 'keep_odd']
    if pairing not in valid_pairings:
        raise ValueError('pairing must be one of %s, not %s'
                         % (valid_pairings, pairing))
    if len(z) == len(p) == 0:
        return array([[k, 0., 0., 1., 0., 0.]])

    # ensure we have the same number of poles and zeros, and make copies
    p = np.concatenate((p, np.zeros(max(len(z) - len(p), 0))))
    z = np.concatenate((z, np.zeros(max(len(p) - len(z), 0))))
    n_sections = (max(len(p), len(z)) + 1) // 2
    sos = zeros((n_sections, 6))

    if len(p) % 2 == 1 and pairing == 'nearest':
        p = np.concatenate((p, [0.]))
        z = np.concatenate((z, [0.]))
    assert len(p) == len(z)

    # Ensure we have complex conjugate pairs
    # (note that _cplxreal only gives us one element of each complex pair):
    z = np.concatenate(_cplxreal(z))
    p = np.concatenate(_cplxreal(p))

    p_sos = np.zeros((n_sections, 2), np.complex128)
    z_sos = np.zeros_like(p_sos)
    for si in range(n_sections):
        # Select the next "worst" pole
        p1_idx = np.argmin(np.abs(1 - np.abs(p)))
        p1 = p[p1_idx]
        p = np.delete(p, p1_idx)

        # Pair that pole with a zero

        if np.isreal(p1) and np.isreal(p).sum() == 0:
            # Special case to set a first-order section
            z1_idx = _nearest_real_complex_idx(z, p1, 'real')
            z1 = z[z1_idx]
            z = np.delete(z, z1_idx)
            p2 = z2 = 0
        else:
            if not np.isreal(p1) and np.isreal(z).sum() == 1:
                # Special case to ensure we choose a complex zero to pair
                # with so later (setting up a first-order section)
                z1_idx = _nearest_real_complex_idx(z, p1, 'complex')
                assert not np.isreal(z[z1_idx])
            else:
                # Pair the pole with the closest zero (real or complex)
                z1_idx = np.argmin(np.abs(p1 - z))
            z1 = z[z1_idx]
            z = np.delete(z, z1_idx)

            # Now that we have p1 and z1, figure out what p2 and z2 need to be
            if not np.isreal(p1):
                if not np.isreal(z1):  # complex pole, complex zero
                    p2 = p1.conj()
                    z2 = z1.conj()
                else:  # complex pole, real zero
                    p2 = p1.conj()
                    z2_idx = _nearest_real_complex_idx(z, p1, 'real')
                    z2 = z[z2_idx]
                    assert np.isreal(z2)
                    z = np.delete(z, z2_idx)
            else:
                if not np.isreal(z1):  # real pole, complex zero
                    z2 = z1.conj()
                    p2_idx = _nearest_real_complex_idx(p, z1, 'real')
                    p2 = p[p2_idx]
                    assert np.isreal(p2)
                else:  # real pole, real zero
                    # pick the next "worst" pole to use
                    idx = np.where(np.isreal(p))[0]
                    assert len(idx) > 0
                    p2_idx = idx[np.argmin(np.abs(np.abs(p[idx]) - 1))]
                    p2 = p[p2_idx]
                    # find a real zero to match the added pole
                    assert np.isreal(p2)
                    z2_idx = _nearest_real_complex_idx(z, p2, 'real')
                    z2 = z[z2_idx]
                    assert np.isreal(z2)
                    z = np.delete(z, z2_idx)
                p = np.delete(p, p2_idx)
        p_sos[si] = [p1, p2]
        z_sos[si] = [z1, z2]
    assert len(p) == len(z) == 0  # we've consumed all poles and zeros
    del p, z

    # Construct the system, reversing order so the "worst" are last
    p_sos = np.reshape(p_sos[::-1], (n_sections, 2))
    z_sos = np.reshape(z_sos[::-1], (n_sections, 2))
    gains = np.ones(n_sections)
    gains[0] = k
    for si in range(n_sections):
        x = zpk2tf(z_sos[si], p_sos[si], gains[si])
        sos[si] = np.concatenate(x)
    return sos

def sosfilt(sos, x, axis=-1, zi=None):
    """Stolen from scipy, please kill me and upgrade scipy...
    """

    import numpy as np
    from numpy import atleast_1d, atleast_2d, array
    from scipy.signal import lfilter

    x = np.asarray(x)

    sos = atleast_2d(sos)
    if sos.ndim != 2:
        raise ValueError('sos array must be 2D')

    n_sections, m = sos.shape
    if m != 6:
        raise ValueError('sos array must be shape (n_sections, 6)')

    use_zi = zi is not None
    if use_zi:
        zi = np.asarray(zi)
        e_zi_shape = list(x.shape)
        x_zi_shape[axis] = 2
        x_zi_shape = tuple([n_sections] + x_zi_shape)
        if zi.shape != x_zi_shape:
            raise ValueError('Invalid zi shape.  With axis=%r, an input with '
                             'shape %r, and an sos array with %d sections, zi '
                             'must have shape %r.' %
                             (axis, x.shape, n_sections, x_zi_shape))
        zf = zeros_like(zi)

    for section in range(n_sections):
        if use_zi:
            x, zf[section] = lfilter(sos[section, :3], sos[section, 3:],
                                     x, axis, zi=zi[section])
        else:
            x = lfilter(sos[section, :3], sos[section, 3:], x, axis)
    out = (x, zf) if use_zi else x
    return out

def _cplxreal(z, tol=None):

    import numpy as np
    from numpy import atleast_1d, atleast_2d, array

    z = atleast_1d(z)
    if z.size == 0:
        return z, z
    elif z.ndim != 1:
        raise ValueError('_cplxreal only accepts 1D input')

    if tol is None:
        # Get tolerance from dtype of input
        tol = 100 * np.finfo((1.0 * z).dtype).eps

    # Sort by real part, magnitude of imaginary part (speed up further sorting)
    z = z[np.lexsort((abs(z.imag), z.real))]

    # Split reals from conjugate pairs
    real_indices = abs(z.imag) <= tol * abs(z)
    zr = z[real_indices].real

    if len(zr) == len(z):
        # Input is entirely real
        return array([]), zr

    # Split positive and negative halves of conjugates
    z = z[~real_indices]
    zp = z[z.imag > 0]
    zn = z[z.imag < 0]

    if len(zp) != len(zn):
        raise ValueError('Array contains complex value with no matching '
                         'conjugate.')

    # Find runs of (approximately) the same real part
    same_real = np.diff(zp.real) <= tol * abs(zp[:-1])
    diffs = numpy.diff(concatenate(([0], same_real, [0])))
    run_starts = numpy.where(diffs > 0)[0]
    run_stops = numpy.where(diffs < 0)[0]

    # Sort each run by their imaginary parts
    for i in range(len(run_starts)):
        start = run_starts[i]
        stop = run_stops[i] + 1
        for chunk in (zp[start:stop], zn[start:stop]):
            chunk[...] = chunk[np.lexsort([abs(chunk.imag)])]

    # Check that negatives match positives
    if any(abs(zp - zn.conj()) > tol * abs(zn)):
        raise ValueError('Array contains complex value with no matching '
                         'conjugate.')

    # Average out numerical inaccuracy in real vs imag parts of pairs
    zc = (zp + zn.conj()) / 2

    return zc, zr

def _nearest_real_complex_idx(fro, to, which):

    import numpy as np

    """Get the next closest real or complex element based on distance"""
    assert which in ('real', 'complex')
    order = np.argsort(np.abs(fro - to))
    mask = np.isreal(fro[order])
    if which == 'complex':
        mask = ~mask
    return order[np.where(mask)[0][0]]

"""A parser for HTML and XHTML."""

# This file is based on sgmllib.py, but the API is slightly different.

# XXX There should be a way to distinguish between PCDATA (parsed
# character data -- the normal case), RCDATA (replaceable character
# data -- only char and entity references and end tags are special)
# and CDATA (character data -- only end tags are special).


import markupbase
import re

# Regular expressions used for parsing

interesting_normal = re.compile('[&<]')
incomplete = re.compile('&[a-zA-Z#]')

entityref = re.compile('&([a-zA-Z][-.a-zA-Z0-9]*)[^a-zA-Z0-9]')
charref = re.compile('&#(?:[0-9]+|[xX][0-9a-fA-F]+)[^0-9a-fA-F]')

starttagopen = re.compile('<[a-zA-Z]')
piclose = re.compile('>')
commentclose = re.compile(r'--\s*>')

# see http://www.w3.org/TR/html5/tokenization.html#tag-open-state
# and http://www.w3.org/TR/html5/tokenization.html#tag-name-state
# note: if you change tagfind/attrfind remember to update locatestarttagend too
tagfind = re.compile('([a-zA-Z][^\t\n\r\f />\x00]*)(?:\s|/(?!>))*')
# this regex is currently unused, but left for backward compatibility
tagfind_tolerant = re.compile('[a-zA-Z][^\t\n\r\f />\x00]*')

attrfind = re.compile(
    r'((?<=[\'"\s/])[^\s/>][^\s/=>]*)(\s*=+\s*'
    r'(\'[^\']*\'|"[^"]*"|(?![\'"])[^>\s]*))?(?:\s|/(?!>))*')

locatestarttagend = re.compile(r"""
  <[a-zA-Z][^\t\n\r\f />\x00]*       # tag name
  (?:[\s/]*                          # optional whitespace before attribute name
    (?:(?<=['"\s/])[^\s/>][^\s/=>]*  # attribute name
      (?:\s*=+\s*                    # value indicator
        (?:'[^']*'                   # LITA-enclosed value
          |"[^"]*"                   # LIT-enclosed value
          |(?!['"])[^>\s]*           # bare value
         )
       )?(?:\s|/(?!>))*
     )*
   )?
  \s*                                # trailing whitespace
""", re.VERBOSE)
endendtag = re.compile('>')
# the HTML 5 spec, section 8.1.2.2, doesn't allow spaces between
# </ and the tag name, so maybe this should be fixed
endtagfind = re.compile('</\s*([a-zA-Z][-.a-zA-Z0-9:_]*)\s*>')


class HTMLParseError(Exception):
    """Exception raised for all parse errors."""

    def __init__(self, msg, position=(None, None)):
        assert msg
        self.msg = msg
        self.lineno = position[0]
        self.offset = position[1]

    def __str__(self):
        result = self.msg
        if self.lineno is not None:
            result = result + ", at line %d" % self.lineno
        if self.offset is not None:
            result = result + ", column %d" % (self.offset + 1)
        return result


class HTMLParser(markupbase.ParserBase):
    """Find tags and other markup and call handler functions.

    Usage:
        p = HTMLParser()
        p.feed(data)
        ...
        p.close()

    Start tags are handled by calling self.handle_starttag() or
    self.handle_startendtag(); end tags by self.handle_endtag().  The
    data between tags is passed from the parser to the derived class
    by calling self.handle_data() with the data as argument (the data
    may be split up in arbitrary chunks).  Entity references are
    passed by calling self.handle_entityref() with the entity
    reference as the argument.  Numeric character references are
    passed to self.handle_charref() with the string containing the
    reference as the argument.
    """

    CDATA_CONTENT_ELEMENTS = ("script", "style")


    def __init__(self):
        """Initialize and reset this instance."""
        self.reset()

    def reset(self):
        """Reset this instance.  Loses all unprocessed data."""
        self.rawdata = ''
        self.lasttag = '???'
        self.interesting = interesting_normal
        self.cdata_elem = None
        markupbase.ParserBase.reset(self)

    def feed(self, data):
        r"""Feed data to the parser.

        Call this as often as you want, with as little or as much text
        as you want (may include '\n').
        """
        self.rawdata = self.rawdata + data
        self.goahead(0)

    def close(self):
        """Handle any buffered data."""
        self.goahead(1)

    def error(self, message):
        raise HTMLParseError(message, self.getpos())

    __starttag_text = None

    def get_starttag_text(self):
        """Return full source of start tag: '<...>'."""
        return self.__starttag_text

    def set_cdata_mode(self, elem):
        self.cdata_elem = elem.lower()
        self.interesting = re.compile(r'</\s*%s\s*>' % self.cdata_elem, re.I)

    def clear_cdata_mode(self):
        self.interesting = interesting_normal
        self.cdata_elem = None

    # Internal -- handle data as far as reasonable.  May leave state
    # and data to be processed by a subsequent call.  If 'end' is
    # true, force handling all data as if followed by EOF marker.
    def goahead(self, end):
        rawdata = self.rawdata
        i = 0
        n = len(rawdata)
        while i < n:
            match = self.interesting.search(rawdata, i) # < or &
            if match:
                j = match.start()
            else:
                if self.cdata_elem:
                    break
                j = n
            if i < j: self.handle_data(rawdata[i:j])
            i = self.updatepos(i, j)
            if i == n: break
            startswith = rawdata.startswith
            if startswith('<', i):
                if starttagopen.match(rawdata, i): # < + letter
                    k = self.parse_starttag(i)
                elif startswith("</", i):
                    k = self.parse_endtag(i)
                elif startswith("<!--", i):
                    k = self.parse_comment(i)
                elif startswith("<?", i):
                    k = self.parse_pi(i)
                elif startswith("<!", i):
                    k = self.parse_html_declaration(i)
                elif (i + 1) < n:
                    self.handle_data("<")
                    k = i + 1
                else:
                    break
                if k < 0:
                    if not end:
                        break
                    k = rawdata.find('>', i + 1)
                    if k < 0:
                        k = rawdata.find('<', i + 1)
                        if k < 0:
                            k = i + 1
                    else:
                        k += 1
                    self.handle_data(rawdata[i:k])
                i = self.updatepos(i, k)
            elif startswith("&#", i):
                match = charref.match(rawdata, i)
                if match:
                    name = match.group()[2:-1]
                    self.handle_charref(name)
                    k = match.end()
                    if not startswith(';', k-1):
                        k = k - 1
                    i = self.updatepos(i, k)
                    continue
                else:
                    if ";" in rawdata[i:]:  # bail by consuming '&#'
                        self.handle_data(rawdata[i:i+2])
                        i = self.updatepos(i, i+2)
                    break
            elif startswith('&', i):
                match = entityref.match(rawdata, i)
                if match:
                    name = match.group(1)
                    self.handle_entityref(name)
                    k = match.end()
                    if not startswith(';', k-1):
                        k = k - 1
                    i = self.updatepos(i, k)
                    continue
                match = incomplete.match(rawdata, i)
                if match:
                    # match.group() will contain at least 2 chars
                    if end and match.group() == rawdata[i:]:
                        self.error("EOF in middle of entity or char ref")
                    # incomplete
                    break
                elif (i + 1) < n:
                    # not the end of the buffer, and can't be confused
                    # with some other construct
                    self.handle_data("&")
                    i = self.updatepos(i, i + 1)
                else:
                    break
            else:
                assert 0, "interesting.search() lied"
        # end while
        if end and i < n and not self.cdata_elem:
            self.handle_data(rawdata[i:n])
            i = self.updatepos(i, n)
        self.rawdata = rawdata[i:]

    # Internal -- parse html declarations, return length or -1 if not terminated
    # See w3.org/TR/html5/tokenization.html#markup-declaration-open-state
    # See also parse_declaration in _markupbase
    def parse_html_declaration(self, i):
        rawdata = self.rawdata
        if rawdata[i:i+2] != '<!':
            self.error('unexpected call to parse_html_declaration()')
        if rawdata[i:i+4] == '<!--':
            # this case is actually already handled in goahead()
            return self.parse_comment(i)
        elif rawdata[i:i+3] == '<![':
            return self.parse_marked_section(i)
        elif rawdata[i:i+9].lower() == '<!doctype':
            # find the closing >
            gtpos = rawdata.find('>', i+9)
            if gtpos == -1:
                return -1
            self.handle_decl(rawdata[i+2:gtpos])
            return gtpos+1
        else:
            return self.parse_bogus_comment(i)

    # Internal -- parse bogus comment, return length or -1 if not terminated
    # see http://www.w3.org/TR/html5/tokenization.html#bogus-comment-state
    def parse_bogus_comment(self, i, report=1):
        rawdata = self.rawdata
        if rawdata[i:i+2] not in ('<!', '</'):
            self.error('unexpected call to parse_comment()')
        pos = rawdata.find('>', i+2)
        if pos == -1:
            return -1
        if report:
            self.handle_comment(rawdata[i+2:pos])
        return pos + 1

    # Internal -- parse processing instr, return end or -1 if not terminated
    def parse_pi(self, i):
        rawdata = self.rawdata
        assert rawdata[i:i+2] == '<?', 'unexpected call to parse_pi()'
        match = piclose.search(rawdata, i+2) # >
        if not match:
            return -1
        j = match.start()
        self.handle_pi(rawdata[i+2: j])
        j = match.end()
        return j

    # Internal -- handle starttag, return end or -1 if not terminated
    def parse_starttag(self, i):
        self.__starttag_text = None
        endpos = self.check_for_whole_start_tag(i)
        if endpos < 0:
            return endpos
        rawdata = self.rawdata
        self.__starttag_text = rawdata[i:endpos]

        # Now parse the data between i+1 and j into a tag and attrs
        attrs = []
        match = tagfind.match(rawdata, i+1)
        assert match, 'unexpected call to parse_starttag()'
        k = match.end()
        self.lasttag = tag = match.group(1).lower()

        while k < endpos:
            m = attrfind.match(rawdata, k)
            if not m:
                break
            attrname, rest, attrvalue = m.group(1, 2, 3)
            if not rest:
                attrvalue = None
            elif attrvalue[:1] == '\'' == attrvalue[-1:] or \
                 attrvalue[:1] == '"' == attrvalue[-1:]:
                attrvalue = attrvalue[1:-1]
            if attrvalue:
                attrvalue = self.unescape(attrvalue)
            attrs.append((attrname.lower(), attrvalue))
            k = m.end()

        end = rawdata[k:endpos].strip()
        if end not in (">", "/>"):
            lineno, offset = self.getpos()
            if "\n" in self.__starttag_text:
                lineno = lineno + self.__starttag_text.count("\n")
                offset = len(self.__starttag_text) \
                         - self.__starttag_text.rfind("\n")
            else:
                offset = offset + len(self.__starttag_text)
            self.handle_data(rawdata[i:endpos])
            return endpos
        if end.endswith('/>'):
            # XHTML-style empty tag: <span attr="value" />
            self.handle_startendtag(tag, attrs)
        else:
            self.handle_starttag(tag, attrs)
            if tag in self.CDATA_CONTENT_ELEMENTS:
                self.set_cdata_mode(tag)
        return endpos

    # Internal -- check to see if we have a complete starttag; return end
    # or -1 if incomplete.
    def check_for_whole_start_tag(self, i):
        rawdata = self.rawdata
        m = locatestarttagend.match(rawdata, i)
        if m:
            j = m.end()
            next = rawdata[j:j+1]
            if next == ">":
                return j + 1
            if next == "/":
                if rawdata.startswith("/>", j):
                    return j + 2
                if rawdata.startswith("/", j):
                    # buffer boundary
                    return -1
                # else bogus input
                self.updatepos(i, j + 1)
                self.error("malformed empty start tag")
            if next == "":
                # end of input
                return -1
            if next in ("abcdefghijklmnopqrstuvwxyz=/"
                        "ABCDEFGHIJKLMNOPQRSTUVWXYZ"):
                # end of input in or before attribute value, or we have the
                # '/' from a '/>' ending
                return -1
            if j > i:
                return j
            else:
                return i + 1
        raise AssertionError("we should not get here!")

    # Internal -- parse endtag, return end or -1 if incomplete
    def parse_endtag(self, i):
        rawdata = self.rawdata
        assert rawdata[i:i+2] == "</", "unexpected call to parse_endtag"
        match = endendtag.search(rawdata, i+1) # >
        if not match:
            return -1
        gtpos = match.end()
        match = endtagfind.match(rawdata, i) # </ + tag + >
        if not match:
            if self.cdata_elem is not None:
                self.handle_data(rawdata[i:gtpos])
                return gtpos
            # find the name: w3.org/TR/html5/tokenization.html#tag-name-state
            namematch = tagfind.match(rawdata, i+2)
            if not namematch:
                # w3.org/TR/html5/tokenization.html#end-tag-open-state
                if rawdata[i:i+3] == '</>':
                    return i+3
                else:
                    return self.parse_bogus_comment(i)
            tagname = namematch.group(1).lower()
            # consume and ignore other stuff between the name and the >
            # Note: this is not 100% correct, since we might have things like
            # </tag attr=">">, but looking for > after tha name should cover
            # most of the cases and is much simpler
            gtpos = rawdata.find('>', namematch.end())
            self.handle_endtag(tagname)
            return gtpos+1

        elem = match.group(1).lower() # script or style
        if self.cdata_elem is not None:
            if elem != self.cdata_elem:
                self.handle_data(rawdata[i:gtpos])
                return gtpos

        self.handle_endtag(elem)
        self.clear_cdata_mode()
        return gtpos

    # Overridable -- finish processing of start+end tag: <tag.../>
    def handle_startendtag(self, tag, attrs):
        self.handle_starttag(tag, attrs)
        self.handle_endtag(tag)

    # Overridable -- handle start tag
    def handle_starttag(self, tag, attrs):
        pass

    # Overridable -- handle end tag
    def handle_endtag(self, tag):
        pass

    # Overridable -- handle character reference
    def handle_charref(self, name):
        pass

    # Overridable -- handle entity reference
    def handle_entityref(self, name):
        pass

    # Overridable -- handle data
    def handle_data(self, data):
        pass

    # Overridable -- handle comment
    def handle_comment(self, data):
        pass

    # Overridable -- handle declaration
    def handle_decl(self, decl):
        pass

    # Overridable -- handle processing instruction
    def handle_pi(self, data):
        pass

    def unknown_decl(self, data):
        pass

    # Internal -- helper to remove special character quoting
    entitydefs = None
    def unescape(self, s):
        if '&' not in s:
            return s
        def replaceEntities(s):
            s = s.groups()[0]
            try:
                if s[0] == "#":
                    s = s[1:]
                    if s[0] in ['x','X']:
                        c = int(s[1:], 16)
                    else:
                        c = int(s)
                    return unichr(c)
            except ValueError:
                return '&#'+s+';'
            else:
                # Cannot use name2codepoint directly, because HTMLParser supports apos,
                # which is not part of HTML 4
                import htmlentitydefs
                if HTMLParser.entitydefs is None:
                    entitydefs = HTMLParser.entitydefs = {'apos':u"'"}
                    for k, v in htmlentitydefs.name2codepoint.iteritems():
                        entitydefs[k] = unichr(v)
                try:
                    return self.entitydefs[s]
                except KeyError:
                    return '&'+s+';'

        return re.sub(r"&(#?[xX]?(?:[0-9a-fA-F]+|\w{1,8}));", replaceEntities, s)
