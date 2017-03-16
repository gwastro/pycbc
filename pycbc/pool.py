""" Tools for creating pools of worker processes
"""
from __future__ import absolute_import
import multiprocessing.pool
import types


_process_lock = None    
_numdone = None    
def _lockstep_fcn(values):
    """ Wrapper to ensure that all processes execute together """
    global _numdone
    numrequired, fcn, args = values
    with _process_lock:
        _numdone.value += 1
    # yep this is an ugly busy loop, do something better please
    # when we care about the performance of this call and not just the 
    # guarantee it provides (ok... maybe never)
    while 1:
        if _numdone.value == numrequired:
            return fcn(args)

class BroadcastPool(multiprocessing.pool.Pool):
    """ Multiprocessing pool with a broadcast method
    """
    def __init__(*args, **kwds):
        global _process_lock
        global _numdone
        _process_lock = multiprocessing.Lock()
        _numdone = multiprocessing.Value('i', 0)
        multiprocessing.pool.Pool.__init__(*args, **kwds)

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
        global _numdone
        results = self.map(_lockstep_fcn, [(len(self), fcn, args)] * len(self))
        _numdone.value = 0
        return results

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


