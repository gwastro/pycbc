# Colored log, requires Python 2.3 or up.

import os
import sys
from distutils.log import *
from distutils.log import Log as old_Log
from distutils.log import _global_log

# Hooks for colored terminal output.
# See also http://www.livinglogic.de/Python/ansistyle
def terminal_has_colors():
    if sys.platform=='cygwin' and 'USE_COLOR' not in os.environ:
        # Avoid importing curses that causes illegal operation
        # with a message:
        #  PYTHON2 caused an invalid page fault in
        #  module CYGNURSES7.DLL as 015f:18bbfc28
        # Details: Python 2.3.3 [GCC 3.3.1 (cygming special)]
        #          ssh to Win32 machine from debian
        #          curses.version is 2.2
        #          CYGWIN_98-4.10, release 1.5.7(0.109/3/2))
        return 0
    if hasattr(sys.stdout,'isatty') and sys.stdout.isatty():
        try:
            import curses
            curses.setupterm()
            if (curses.tigetnum("colors") >= 0
                and curses.tigetnum("pairs") >= 0
                and ((curses.tigetstr("setf") is not None
                      and curses.tigetstr("setb") is not None)
                     or (curses.tigetstr("setaf") is not None
                         and curses.tigetstr("setab") is not None)
                     or curses.tigetstr("scp") is not None)):
                return 1
        except Exception,msg:
            pass
    return 0

if terminal_has_colors():
    _colour_codes = dict(black=0, red=1, green=2, yellow=3,
                         blue=4, magenta=5, cyan=6, white=7, default=9)
    def colour_text(s, fg=None, bg=None, bold=False):
        seq = []
        if bold:
            seq.append('1')
        if fg:
            fgcode = 30 + _colour_codes.get(fg.lower(), 0)
            seq.append(str(fgcode))
        if bg:
            bgcode = 40 + _colour_codes.get(fg.lower(), 7)
            seq.append(str(bgcode))
        if seq:
            return '\x1b[%sm%s\x1b[0m' % (';'.join(seq), s)
        else:
            return s
else:
    def colour_text(s, fg=None, bg=None):
        return s

def default_text(s):
    return colour_text(s, 'default')
def red_text(s):
    return colour_text(s, 'red')
def green_text(s):
    return colour_text(s, 'green')
def yellow_text(s):
    return colour_text(s, 'yellow')
def cyan_text(s):
    return colour_text(s, 'cyan')
def blue_text(s):
    return colour_text(s, 'blue')


def is_string(s):
    return isinstance(s, str)

def is_sequence(seq):
    if is_string(seq):
        return False
    try:
        len(seq)
    except:
        return False
    return True


def _fix_args(args,flag=1):
    if is_string(args):
        return args.replace('%','%%')
    if flag and is_sequence(args):
        return tuple([_fix_args(a,flag=0) for a in args])
    return args

class Log(old_Log):
    def _log(self, level, msg, args):
        if level >= self.threshold:
            if args:
                msg = msg % _fix_args(args)
            if 0:
                if msg.startswith('copying ') and msg.find(' -> ') != -1:
                    return
                if msg.startswith('byte-compiling '):
                    return
            print _global_color_map[level](msg)
            sys.stdout.flush()

    def good(self, msg, *args):
        """If we'd log WARN messages, log this message as a 'nice' anti-warn
        message.
        """
        if WARN >= self.threshold:
            if args:
                print green_text(msg % _fix_args(args))
            else:
                print green_text(msg)
            sys.stdout.flush()
_global_log.__class__ = Log

good = _global_log.good

def set_threshold(level, force=False):
    prev_level = _global_log.threshold
    if prev_level > DEBUG or force:
        # If we're running at DEBUG, don't change the threshold, as there's
        # likely a good reason why we're running at this level.
        _global_log.threshold = level
        if level <= DEBUG:
            info('set_threshold: setting thershold to DEBUG level, it can be changed only with force argument')
    else:
        info('set_threshold: not changing thershold from DEBUG level %s to %s' % (prev_level,level))
    return prev_level

def set_verbosity(v, force=False):
    prev_level = _global_log.threshold
    if v < 0:
        set_threshold(ERROR, force)
    elif v == 0:
        set_threshold(WARN, force)
    elif v == 1:
        set_threshold(INFO, force)
    elif v >= 2:
        set_threshold(DEBUG, force)
    return {FATAL:-2,ERROR:-1,WARN:0,INFO:1,DEBUG:2}.get(prev_level,1)

_global_color_map = {
    DEBUG:cyan_text,
    INFO:yellow_text,
    WARN:red_text,
    ERROR:red_text,
    FATAL:red_text
}

# don't use INFO,.. flags in set_verbosity, these flags are for set_threshold.
set_verbosity(0, force=True)
