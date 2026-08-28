
### marg-reconstruct-cache: multiprocessing attempt, not working (2026-08-27)

Alex wants multiprocessing covered even if MPI is not. Attempted via
`BroadcastPool.broadcast`, which runs a function on every worker — a
`gather_stats_cache(pool)` method that merges each worker's store into the
parent after the run. **Reverted, it does not work:** the gather runs and
collects 148 entries where it previously collected none, but the hit rate on
the posterior stays at 0%.

Two things wrong, and they are separate:

- **148 is far too few.** Two workers between them made thousands of
  evaluations with a 200000-entry limit, so most of what they computed is not
  in the stores being returned. That points at the workers not holding the
  model instance the gather reaches — likely `models._global_instance` in a
  worker not being the object whose `update` ran, or the pool having forked
  before it was set.
- **0% hits on 148 real entries** means the keys do not match either, so even
  the gathered stats are filed under something the parent does not ask for.
  Worth checking whether a worker's `update` receives the same parameters the
  parent later looks up, or transformed ones.

Both are worth an hour with a print statement in a worker; neither is a
design flaw in the cache itself. The pushed branch `b9e43f1b5` remains the
working in-process version, +59/-1, and is unaffected.
