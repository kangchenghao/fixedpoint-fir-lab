# Fixed-point FIR Lab

**Project led by [kangchenghaokang-dot](https://github.com/kangchenghaokang-dot).** I directed the wider
research effort toward practical, lightweight and reproducible public work.
This artifact makes that goal concrete through an explicit numerical contract,
comparisons and runnable evidence. See [my contributions and the implementation
division](CONTRIBUTIONS.md).

A small numerical experiment about a practical question: **how much error comes from input quantization, coefficient quantization, and the final output shift?** It uses Python 3.10+ and the standard library. All signals are synthetic. There is no EEG recording, trained model, clinical claim, hardware implementation, or new filtering algorithm.

The experiment preserves failures as well as improvements. In the supplied small-amplitude example, the coarse output fails an illustrative 1% error line; finer output scales reduce its error. That is a numerical result for these constructed inputs, not a universal precision recommendation.

## Run without installing anything

From this repository directory:

```text
python tools/check.py
```

This runs all tests and regenerates `results/demo.json` and `results/synthetic-fixtures.json`. It sets `PYTHONPATH` only for its subprocesses and does not install dependencies or change your system settings. To run only the experiment with no environment setup:

```text
python run.py experiment --output results/demo.json
```

For the requested module CLI, expose the `src` directory first:

```powershell
# PowerShell, in this repository
$env:PYTHONPATH = (Resolve-Path src).Path
python -m unittest discover -s tests -v
python -m fixedpoint_fir_lab experiment --output results/demo.json
```

```sh
# POSIX shell
PYTHONPATH=src python -m fixedpoint_fir_lab experiment --output results/demo.json
```

An installed package also supports the same module command and the `fixedpoint-fir-lab` entry point. Packaging metadata is provided in `pyproject.toml`; installation is optional and was not needed for the checked-in experiment. Building a package uses setuptools; numerical execution has no third-party dependency.

## Numerical contract

The filter is the standard 31-tap, symmetric Hamming-windowed difference of two low-pass sinc functions, with sample rate 160 Hz and design cutoffs 8 and 30 Hz:

```text
m = k - 15
sinc(x) = sin(pi*x)/(pi*x), with sinc(0) = 1
h_raw[k] = (2*30/160*sinc(2*30/160*m)
          - 2*8/160*sinc(2*8/160*m))
           * (0.54 - 0.46*cos(2*pi*k/30))
```

Normalize its gain to one at the passband center, 19 Hz. The cutoff arguments are transition centers (approximately −6 dB here), not brick-wall boundaries. No SciPy routine is called. The [standard window-method documentation](https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.firwin.html) is a useful reference for these conventions; this repository does not claim an external-library conformance test.

Input integers represent quarter-microvolt units: `x_int = round(x_uv * 4)`. Coefficients are `h_q15 = round(h * 32768)`. Both use ties to even. Inputs outside the signed 14-bit domain `[-8192, 8191]` after rounding are rejected rather than silently clipped. Q15 coefficients must fit signed 16-bit integers.

Each exact accumulator is the causal convolution of these integers. All samples before and after the finite input are zero. **N input samples produce N+30 outputs**, retaining startup and tail; no warmup samples are discarded from any metric.

```text
y_int[n] = accumulator[n] >> shift
units_per_uv = 4 * 32768 / 2**shift
y_uv[n] = y_int[n] / units_per_uv
```

| Right shift | Integer units per uV | Output LSB |
|---|---:|---:|
| 15 | 4 | 0.25 uV |
| 13 | 16 | 0.0625 uV |
| 11 | 64 | 0.015625 uV |

Right shift floors toward negative infinity, including negative accumulators; it does not round to nearest or truncate toward zero. This introduces a negative output-stage bias, which the JSON reports. Output integers at different scales are never compared directly.

The linear-phase group delay is 15 samples, or 93.75 ms at 160 Hz. **This is signal-processing group delay, not hardware computation latency.** No hardware timing is measured.

## What is checked

The direct form and symmetric form produce exactly equal integer accumulators and outputs. The symmetric form combines mirrored input samples before multiplication; it requires exactly symmetric integer coefficients. An additional input-major polynomial-product reference validates supplied output and rejects deliberately corrupted values or missing tail samples.

For every direct or paired partial sum, the triangle inequality gives:

```text
abs(partial_sum) <= max_abs_input * sum(abs(h_q15))
                 <= 8192 * 47541 = 389455872 < 2**31
```

The main arithmetic uses arbitrary-precision Python integers. The experiment additionally checks each partial sum with a signed 32-bit overflow guard. Tests cover positive and negative overflow using deliberately unsafe coefficient sets; there is no silent wraparound or saturation. This bound is about accumulator range, not the cost or latency of a fixed-width implementation.

The 23 unittest methods cover independent hand calculations, impulse response, causality, complete tail, frequency-response conventions, linear phase, positive and negative tie cases, input rejection, negative-shift semantics, extreme and seeded inputs, direct/symmetric agreement, overflow, output corruption, metrics, hashes and the actual CLI. These are ordinary executable tests and a simple mathematical bound, not a mechanized formal proof.

## Measured example

Checked-in `results/demo.json` was generated with Python 3.13.6, seed 20261112, and two 512-sample inputs. There are 542 outputs per record. The smooth record combines in-band and out-of-band sinusoids with small seeded noise. The smaller record is deliberately generated on a quarter-uV grid to isolate output precision; its input-quantization error is therefore zero by construction. It is not a held-out biological dataset.

Total NRMSE below is relative to the full convolution of **unquantized design taps and the original synthetic samples**, with all results expressed in uV:

| Record | Shift 15 | Shift 13 | Shift 11 |
|---|---:|---:|---:|
| smooth | 0.8665% | 0.2785% | 0.2039% |
| small_amplitude | 3.8337% | 0.9662% | 0.2423% |

The 1% line is illustrative, not a prevalidated application acceptance criterion. Do not infer that shift 13 works for all signals. The finer output does not remove input or coefficient error.

The JSON separately compares:

1. Quantized input + original taps against the original floating design output.
2. Quantized input + quantized taps against stage 1.
3. Shifted integer output, converted to uV, against stage 2.
4. Final output against the original floating design output.

RMSE values from these stages should not be added or combined as if errors were independent. NRMSE is `100*RMSE/reference_RMS`; it is `null` when reference RMS is zero. The report also includes maximum absolute error, bias, frequency response, the accumulator bound, parameters, fixture hashes, executing package-source hashes and environment metadata. It contains no machine-specific absolute paths.

## Explore further

- Change amplitudes while keeping the filter and scales fixed. Where does output precision stop being the main error source?
- Add a new synthetic test with its purpose declared before inspecting results.
- Compare truncation-to-zero with the current floor rule and inspect bias, especially near zero.
- Vary tap count or cutoffs, then recompute the accumulator bound rather than reusing this filter's number.
- If benchmarking software, report interpreter, method, startup and repetition rules. Algebraic equality alone does not establish a speedup.

See [the short Chinese guide](docs/说明.md), [limitations](docs/LIMITATIONS.md), [sources](SOURCES.md), and [AI assistance](AI_USAGE.md). Own repository code is distributed under the [MIT license](LICENSE); cited documentation and papers keep their own rights.
