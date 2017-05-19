""" Tools for creating pools of worker processes
"""
from __future__ import absolute_import
import multiprocessing.pool
import functools
from multiprocessing import TimeoutError
import types
import signal
import atexit

# Allow the pool to be interrupted, need to disable the children processes
# from intercepting the keyboard interrupt
def _noint(init, *args):
    signal.signal(signal.SIGINT, signal.SIG_IGN)
    if init is not None:
        return init(*args)

_process_lock = None    
_numdone = None    
def _lockstep_fcn(values):
    """ Wrapper to ensure that all processes execute together """
    numrequired, fcn, args = values
    with _process_lock:
        _numdone.value += 1
    # yep this is an ugly busy loop, do something better please
    # when we care about the performance of this call and not just the 
    # guarantee it provides (ok... maybe never)
    while 1:
        if _numdone.value == numrequired:
            return fcn(args)

def _shutdown_pool(p):
    p.terminate()
    p.join()

class BroadcastPool(multiprocessing.pool.Pool):
    """ Multiprocessing pool with a broadcast method
    """
    def __init__(self, processes=None, initializer=None, initargs=(), **kwds):
        global _process_lock
        global _numdone
        _process_lock = multiprocessing.Lock()
        _numdone = multiprocessing.Value('i', 0)
        noint = functools.partial(_noint, initializer)
        super(BroadcastPool, self).__init__(processes, noint, initargs, **kwds)
        atexit.register(_shutdown_pool, self)
        
    def __len__(self):
        return len(self._pool)

    def broadcast(self, fcn, args):
        """ Do a function call on every worker.

        Parameters
        ----------
        fcn: funtion
            Function to call.
        args: tuple
            The arguments for Pool.map
        """
        results = self.map(_lockstep_fcn, [(len(self), fcn, args)] * len(self))
        _numdone.value = 0
        return results

    def map(self, func, items, chunksize=None):
        """ Catch keyboard interuppts to allow the pool to exit cleanly.

        Parameters
        ----------
        func: function
            Function to call
        items: list of tuples
            Arguments to pass
        chunksize: int, Optional
            Number of calls for each process to handle at once
        """
        results = self.map_async(func, items, chunksize)
        while True:
            try:
                return results.get(1800)
            except TimeoutError:
                pass
            except KeyboardInterrupt:
                self.terminate()
                self.join()
                raise KeyboardInterrupt

def _dummy_broadcast(self, f, args):
    self.map(f, [args] * self.size)    

def choose_pool(processes, mpi=False):
    if mpi:
        try:
            import schwimmbad
            pool = schwimmbad.choose_pool(mpi=mpi,
                                          processes=processes)
            pool.broadcast = types.MethodType(_dummy_broadcast, pool)
        except ImportError:
            raise ValueError("Failed to start up an MPI pool, "
                             "install mpi4py / schwimmbadd")
    else:
        pool = BroadcastPool(processes)
    return pool


