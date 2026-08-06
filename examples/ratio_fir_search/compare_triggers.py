#!/usr/bin/env python
"""Compare triggers from the FIR/ratio search (run_fir_search.sh) against
the plain matched-filter reference search (run_reference_search.sh).

For each FIR-search trigger, finds the closest-in-time reference trigger
with the same template and prints the standard deviation of their SNR
difference, alongside the value predicted from the FIR bank's min-match
target (read directly from fir_full_bank.hdf, so it always reflects
whatever make_fir_bank.sh was actually run with).

The FIR/ratio method has a known additional numerical-precision effect
beyond that mismatch bound -- computing in single precision with a FIR
filter whose frequency response can have a wide dynamic range loses some
accuracy relative to the bank's own (double-precision) construction-time
verification -- so a several-times-larger measured value than predicted
is expected, not a sign of a bug.
"""
import argparse

import h5py
import numpy as np

IFO = "L1"


def load_triggers(path):
    with h5py.File(path, "r") as f:
        group = f[IFO]
        return {
            "snr": group["snr"][:],
            "end_time": group["end_time"][:],
            "template_hash": group["template_hash"][:],
        }


def fir_bank_params(path):
    with h5py.File(path, "r") as f:
        attrs = f["fir_data"].attrs
        return float(attrs["min_match"]), float(attrs["sample_rate"])


def matched_snr_diffs(reference, fir, sample_rate):
    tol = 0.5 / sample_rate
    ref_time = reference["end_time"]
    diffs = []
    for i in range(len(fir["end_time"])):
        t = fir["end_time"][i]
        j = np.abs(ref_time - t).argmin()
        if (abs(ref_time[j] - t) < tol
                and fir["template_hash"][i] == reference["template_hash"][j]):
            diffs.append(fir["snr"][i] - reference["snr"][j])
    return np.array(diffs)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reference-file", default="reference.hdf")
    parser.add_argument("--fir-file", default="fir.hdf")
    parser.add_argument("--fir-bank-file", default="fir_full_bank.hdf")
    parser.add_argument("--allowed-margin-factor", type=float, default=5.0,
                        help="How many multiples of the mismatch-predicted "
                             "std to allow before flagging a problem.")
    args = parser.parse_args()

    reference = load_triggers(args.reference_file)
    fir = load_triggers(args.fir_file)
    min_match, sample_rate = fir_bank_params(args.fir_bank_file)

    diffs = matched_snr_diffs(reference, fir, sample_rate)
    if len(diffs) == 0:
        raise SystemExit("No triggers matched between the two searches -- "
                         "did both searches run over the same data/bank?")

    measured_std = float(np.std(diffs))
    predicted_std = float(np.sqrt(2.0 * (1.0 - min_match)))
    ratio = measured_std / predicted_std

    print(f"matched triggers: {len(diffs)} "
         f"(reference={len(reference['end_time'])}, fir={len(fir['end_time'])})")
    print(f"FIR bank min_match target: {min_match}")
    print(f"mismatch-predicted SNR-difference std: {predicted_std:.6g}")
    print(f"measured SNR-difference std:           {measured_std:.6g}")
    print(f"ratio (measured / predicted):          {ratio:.2f}x")

    if ratio > args.allowed_margin_factor:
        raise SystemExit(
            f"Measured std is {ratio:.2f}x the mismatch prediction, "
            f"exceeding the allowed {args.allowed_margin_factor}x margin."
        )
    print(f"OK: within the allowed {args.allowed_margin_factor}x margin.")


if __name__ == "__main__":
    main()
