
## Design decided, not yet built: a generic model-stats cache (2026-08-27)

Branch `marg-reconstruct-cache` (`wt-reccache`, off master `e1904358e`) is
created but **empty** — the first attempt was scrapped because it knew about
marginalization, which is the wrong shape. Alex's steer: *"it should work in
such a way to be clear that it doesn't know about marginalization. It works
just with the extra stats stuff and shouldn't care what is and is not there."*

**What it should be:** a bounded cache of `current_params -> current_stats`
in `BaseModel` (`models/base.py`), which is the class that already owns
`_current_params` and `_current_stats`. It caches whatever the stats hold and
knows nothing about what produced them. A miss costs the evaluation it would
have cost anyway, so correctness never depends on it — which is what makes it
safe under MPI, where every lookup simply misses.

**Why this ordering:** it is the sampler-independent answer to extra stats
being dropped (the reason `pycbc_inference_model_stats` exists), it is
reviewable on its own, and `marg-inline-reconstruct` then reduces to putting
its draws into the stats and getting persistence for free.

**Trap found while attempting it:** a cache key built from `current_params`
must not assume the values are numeric — `approximant = TaylorF2` is in there.

`marg-inline-reconstruct` (`2a1d488e6`, pushed) is unaffected and remains the
reviewable branch for the demarginalization work itself.
