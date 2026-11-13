"""Seeded synthetic experiments and portable evidence, using only stdlib."""

import hashlib
import json
import math
import platform
import random
from pathlib import Path

from .core import (
    accumulator_bound, design_bandpass, fir_accumulators, fir_float, fir_integer,
    frequency_response, quantize_coefficients, quantize_input, verify_integer_output,
)

DEFAULT_SEED = 20261112


def canonical_hash(value):
    data = json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False).encode("utf-8")
    return hashlib.sha256(data).hexdigest()


def metrics(actual, reference):
    if len(actual) != len(reference) or not actual:
        raise ValueError("Metrics require equal nonempty sequences")
    errors = [a - b for a, b in zip(actual, reference)]
    rms = math.sqrt(math.fsum(value * value for value in reference) / len(reference))
    rmse = math.sqrt(math.fsum(value * value for value in errors) / len(errors))
    return {
        "samples": len(errors), "reference_rms_uv": rms,
        "rmse_uv": rmse, "max_abs_error_uv": max(abs(value) for value in errors),
        "bias_uv": math.fsum(errors) / len(errors),
        "nrmse_percent_of_reference_rms": 100 * rmse / rms if rms else None,
    }


def make_fixtures(seed=DEFAULT_SEED, samples=512):
    if type(seed) is not int or not 0 <= seed <= 0xffffffff:
        raise ValueError("seed must be an unsigned 32-bit integer")
    if type(samples) is not int or not 32 <= samples <= 100000:
        raise ValueError("samples must be an integer between 32 and 100000")
    rng = random.Random(seed)
    smooth, small = [], []
    for n in range(samples):
        time = n / 160
        smooth.append(28 * math.sin(2 * math.pi * 12 * time)
                      + 10 * math.sin(2 * math.pi * 22 * time)
                      + 8 * math.sin(2 * math.pi * 2 * time)
                      + 4 * math.sin(2 * math.pi * 55 * time) + rng.uniform(-0.2, 0.2))
        small_raw = (6 * math.sin(2 * math.pi * 12 * time)
                     + 2 * math.sin(2 * math.pi * 22 * time)
                     + math.sin(2 * math.pi * 2 * time) + rng.uniform(-0.025, 0.025))
        # Deliberately grid-aligned synthetic data isolates output precision.
        # It is not an unquantized biological recording or a fitted dataset.
        small.append(round(small_raw * 4) / 4)
    return {"schema": "synthetic-fir-fixtures-v1", "seed": seed, "fs_hz": 160,
            "source": "generated synthetic sinusoids and seeded uniform noise; no human data",
            "records": [
                {"id": "smooth", "grid_policy": "unrounded floating synthetic samples", "samples_uv": smooth},
                {"id": "small_amplitude", "grid_policy": "quarter-uV grid by design to isolate output quantization", "samples_uv": small},
            ]}


def _source_hashes():
    # Installed wheels need not contain the original project directory.
    # Package sources are present in both source-layout and installed runs.
    directory = Path(__file__).resolve().parent
    return {f"src/fixedpoint_fir_lab/{p.name}": hashlib.sha256(p.read_bytes()).hexdigest()
            for p in sorted(directory.glob("*.py"))}


def run_experiment(*, seed=DEFAULT_SEED, samples=512):
    fixtures = make_fixtures(seed, samples)
    taps = design_bandpass()
    q15 = quantize_coefficients(taps)
    quantized_taps = [value / 32768 for value in q15]
    bound = accumulator_bound(q15)
    cases = []
    for record in fixtures["records"]:
        original = record["samples_uv"]
        x = quantize_input(original)
        reference = fir_float(original, taps)
        input_stage = fir_float([value / 4 for value in x], taps)
        coefficient_stage = fir_float([value / 4 for value in x], quantized_taps)
        accumulators = fir_accumulators(x, q15, accumulator_bits=32)
        if coefficient_stage != [value / (4 * 32768) for value in accumulators]:
            raise RuntimeError("Quantized-coefficient floating reference disagrees with exact accumulator")
        variants = []
        for shift in (15, 13, 11):
            output = fir_integer(x, q15, output_right_shift=shift, accumulator_bits=32)
            paired = fir_integer(x, q15, method="symmetric", output_right_shift=shift, accumulator_bits=32)
            if output != paired:
                raise RuntimeError("Direct and symmetric integer results disagree")
            verified = verify_integer_output(x, q15, output, output_right_shift=shift)
            units = 4 * 32768 / (2 ** shift)
            output_uv = [value / units for value in output]
            total = metrics(output_uv, reference)
            relative = total["nrmse_percent_of_reference_rms"]
            variants.append({
                "output_right_shift": shift, "output_units_per_uv": units,
                "output_lsb_uv": 1 / units,
                "output_quantization_incremental": metrics(output_uv, coefficient_stage),
                "total_vs_float_design": total,
                "illustrative_total_nrmse_under_1_percent": relative < 1 if relative is not None else None,
                "output_integer_sha256": canonical_hash(output),
                "direct_symmetric_bitexact": True, "independent_samples_checked": verified,
            })
        cases.append({
            "id": record["id"], "grid_policy": record["grid_policy"],
            "input_samples": len(x), "output_samples": len(reference),
            "samples_uv_sha256": canonical_hash(original), "samples_int_sha256": canonical_hash(x),
            "quantized_input_min": min(x), "quantized_input_max": max(x), "clipped_samples": 0,
            "input_quantization_incremental": metrics(input_stage, reference),
            "coefficient_quantization_incremental": metrics(coefficient_stage, input_stage),
            "input_and_coefficient_total": metrics(coefficient_stage, reference),
            "maximum_observed_abs_accumulator": max(abs(value) for value in accumulators),
            "variants": variants,
        })
    frequencies = [0, 2, 8, 12, 19, 22, 30, 55, 80]
    response = frequency_response(taps, frequencies)
    quantized_response = frequency_response(quantized_taps, frequencies)
    report = {
        "schema": "fixedpoint-fir-experiment-v1", "experiment_complete": True,
        "parameters": {
            "fs_hz": 160, "cutoffs_hz": [8, 30], "taps": 31, "window": "symmetric Hamming",
            "normalization_hz": 19, "input_units_per_uv": 4,
            "input_integer_range": [-8192, 8191], "input_rounding": "ties-to-even; overflow rejected",
            "coefficient_fraction_bits": 15, "coefficient_rounding": "ties-to-even",
            "output_rounding": "arithmetic right shift; floor toward negative infinity",
            "coefficients_float": taps, "coefficients_q15": q15,
            "output_shifts": [15, 13, 11], "seed": seed,
        },
        "accumulator": {"arithmetic": "unbounded Python integers; optional per-partial-sum signed guard",
                        "input_abs_bound": 8192, "sum_abs_coefficients_q15": sum(abs(v) for v in q15),
                        "conservative_abs_bound": bound, "signed_int32_max": 2**31 - 1,
                        "all_direct_and_paired_partial_sums_safe_in_int32": bound <= 2**31 - 1},
        "measurement": {"convolution": "full linear convolution, zero before and after input",
                        "output_length": "N+30", "startup_and_tail_included_in_all_metrics": True,
                        "group_delay_samples": 15, "group_delay_ms": 93.75,
                        "group_delay_is_hardware_compute_latency": False,
                        "threshold_note": "1% is an illustrative numerical comparison, not a clinical or application acceptance criterion",
                        "nrmse_definition": "100 * RMS(actual-reference) / RMS(reference); null when reference RMS is zero",
                        "timing": "no host benchmark, hardware cycles, synthesis, FPGA or power measurements"},
        "frequency_response": [{"frequency_hz": f, "float_magnitude": abs(a), "q15_magnitude": abs(b)}
                               for f, a, b in zip(frequencies, response, quantized_response)],
        "fixtures": {"file": "synthetic-fixtures.json", "canonical_sha256": canonical_hash(fixtures)},
        "cases": cases,
        "sources_sha256": _source_hashes(),
        "environment": {"python": platform.python_version(), "implementation": platform.python_implementation(),
                        "system": platform.system(), "release": platform.release(), "machine": platform.machine()},
        "limitations": ["Synthetic signals only; no EEG recording or clinical validation.",
                        "This is a standard FIR numerical experiment, not a novel filtering algorithm.",
                        "Float design is a numerical reference, not an ideal brick-wall filter.",
                        "One filter and two constructed records do not select universal precision.",
                        "Symmetric arithmetic equality is not a measured speedup.",
                        "Binary floating results may differ slightly across Python/libm platforms; record hashes per run."],
    }
    return report, fixtures


def save_experiment(output, *, seed=DEFAULT_SEED, samples=512):
    target = Path(output)
    if target.suffix.lower() != ".json" or target.name.lower() == "synthetic-fixtures.json":
        raise ValueError("Use a .json report name other than synthetic-fixtures.json")
    report, fixtures = run_experiment(seed=seed, samples=samples)
    target.parent.mkdir(parents=True, exist_ok=True)
    (target.parent / "synthetic-fixtures.json").write_text(
        json.dumps(fixtures, indent=2, allow_nan=False) + "\n", encoding="utf-8", newline="\n")
    target.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n", encoding="utf-8", newline="\n")
    return report
