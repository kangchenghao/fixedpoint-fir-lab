"""Explicit numerical contracts and complete zero-state convolutions.

All convolutions include the startup transient and the zero-input tail, so an
N-sample input and M taps produce N+M-1 outputs. Integer arithmetic is unbounded
unless an optional signed accumulator guard is requested.
"""

import cmath
import math

INPUT_SCALE = 4
INPUT_MIN, INPUT_MAX = -8192, 8191
COEFFICIENT_SCALE = 32768


def _int(value, name, lower, upper=None):
    if type(value) is not int or value < lower or (upper is not None and value > upper):
        raise ValueError(f"{name} must be an integer in [{lower}, {upper}]")


def _real(value, name):
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ValueError(f"{name} must be a finite real number")
    try:
        result = float(value)
    except OverflowError as error:
        raise ValueError(f"{name} is too large") from error
    if not math.isfinite(result):
        raise ValueError(f"{name} must be finite")
    return result


def _sinc(x):
    return 1.0 if x == 0 else math.sin(math.pi * x) / (math.pi * x)


def design_bandpass(fs_hz=160.0, low_hz=8.0, high_hz=30.0, taps=31):
    """Windowed-sinc difference, symmetric Hamming window, center gain = 1.

    This is an independent implementation of the standard design formula.
    Cutoff arguments are design transition centers, not ideal brick-wall edges.
    """
    fs, low, high = (_real(v, n) for v, n in (
        (fs_hz, "fs_hz"), (low_hz, "low_hz"), (high_hz, "high_hz")))
    _int(taps, "taps", 3)
    if taps % 2 != 1 or not (0 < low < high < fs / 2):
        raise ValueError("Require odd taps >= 3 and 0 < low < high < fs/2")
    middle = (taps - 1) // 2
    coefficients = []
    # Calculate one side and mirror exactly, including after normalization.
    for n in range(middle + 1):
        offset = n - middle
        ideal = 2 * high / fs * _sinc(2 * high / fs * offset)
        ideal -= 2 * low / fs * _sinc(2 * low / fs * offset)
        coefficients.append(ideal * (0.54 - 0.46 * math.cos(2 * math.pi * n / (taps - 1))))
    coefficients += coefficients[-2::-1]
    center_hz = (low + high) / 2
    gain = math.fsum(h * math.cos(2 * math.pi * center_hz / fs * (n - middle))
                     for n, h in enumerate(coefficients))
    if abs(gain) < 1e-12:
        raise ValueError("Degenerate normalization gain")
    return [h / gain for h in coefficients]


def quantize_input(samples_uv):
    """uV * 4, ties to even; reject out-of-range values rather than clipping."""
    result = []
    for sample in samples_uv:
        scaled = _real(sample, "sample_uv") * INPUT_SCALE
        if not math.isfinite(scaled):
            raise ValueError("Input is too large")
        value = round(scaled)
        if not INPUT_MIN <= value <= INPUT_MAX:
            raise ValueError("Quantized input exceeds the signed 14-bit domain")
        result.append(value)
    return result


def quantize_coefficients(coefficients):
    """Signed Q15 coefficient integers, using Python round (ties to even)."""
    result = []
    for coefficient in coefficients:
        scaled = _real(coefficient, "coefficient") * COEFFICIENT_SCALE
        if not math.isfinite(scaled):
            raise ValueError("Coefficient is too large")
        value = round(scaled)
        if not -32768 <= value <= 32767:
            raise ValueError("Coefficient exceeds signed Q15 range")
        result.append(value)
    return result


def _inputs(samples, coefficients):
    x, h = list(samples), list(coefficients)
    if not x or not h:
        raise ValueError("Input and coefficient sequences must be nonempty")
    for sample in x:
        _int(sample, "sample_int", INPUT_MIN, INPUT_MAX)
    for coefficient in h:
        _int(coefficient, "coefficient_q15", -32768, 32767)
    return x, h


def accumulator_bound(coefficients, input_bound=8192):
    """Conservative bound for every direct or paired partial sum."""
    _int(input_bound, "input_bound", 0)
    h = list(coefficients)
    if not h:
        raise ValueError("Coefficients must be nonempty")
    for value in h:
        _int(value, "coefficient_q15", -32768, 32767)
    return sum(abs(value) for value in h) * input_bound


def guard_accumulator(value, bits=32):
    """Reject overflow; this is not wrapping or saturating arithmetic."""
    _int(bits, "bits", 2, 256)
    if type(value) is not int or not -(1 << (bits - 1)) <= value < (1 << (bits - 1)):
        raise OverflowError(f"Value exceeds a signed {bits}-bit accumulator")
    return value


def fir_accumulators(samples_int, coefficients_q15, *, method="direct", accumulator_bits=None):
    x, h = _inputs(samples_int, coefficients_q15)
    if method not in ("direct", "symmetric"):
        raise ValueError("method must be direct or symmetric")
    if method == "symmetric" and (len(h) % 2 == 0 or h != h[::-1]):
        raise ValueError("Symmetric mode requires odd, exactly symmetric coefficients")
    if accumulator_bits is not None:
        guard_accumulator(0, accumulator_bits)

    def sample(index):
        return x[index] if 0 <= index < len(x) else 0

    result = []
    for n in range(len(x) + len(h) - 1):
        total = 0
        if method == "direct":
            terms = (coefficient * sample(n - k) for k, coefficient in enumerate(h))
        else:
            half = len(h) // 2
            pairs = [h[k] * (sample(n - k) + sample(n - (len(h) - 1 - k))) for k in range(half)]
            terms = iter(pairs + [h[half] * sample(n - half)])
        for term in terms:
            total += term
            if accumulator_bits is not None:
                guard_accumulator(total, accumulator_bits)
        result.append(total)
    return result


def fir_integer(samples_int, coefficients_q15, *, output_right_shift=13, method="direct", accumulator_bits=None):
    _int(output_right_shift, "output_right_shift", 0, 31)
    return [value >> output_right_shift for value in fir_accumulators(
        samples_int, coefficients_q15, method=method, accumulator_bits=accumulator_bits)]


def fir_float(samples, coefficients):
    x = [_real(v, "sample") for v in samples]
    h = [_real(v, "coefficient") for v in coefficients]
    if not x or not h:
        raise ValueError("Input and coefficient sequences must be nonempty")
    return [math.fsum(h[k] * x[n - k] for k in range(len(h)) if 0 <= n - k < len(x))
            for n in range(len(x) + len(h) - 1)]


def verify_integer_output(samples_int, coefficients_q15, output, *, output_right_shift=13):
    """Independent input-major polynomial product checks supplied output.

    This checks arithmetic only, not timing or independent implementation proof.
    It does not call either production convolution path.
    """
    x, h = _inputs(samples_int, coefficients_q15)
    _int(output_right_shift, "output_right_shift", 0, 31)
    expected = [0] * (len(x) + len(h) - 1)
    for n, sample in enumerate(x):
        for k, coefficient in enumerate(h):
            expected[n + k] += sample * coefficient
    actual = list(output)
    if any(type(v) is not int for v in actual):
        raise ValueError("Output must contain integers")
    expected = [v >> output_right_shift for v in expected]
    if actual != expected:
        raise ValueError("Integer output differs from the input-major reference")
    return len(expected)


def frequency_response(coefficients, frequencies_hz, fs_hz=160.0):
    fs = _real(fs_hz, "fs_hz")
    if fs <= 0:
        raise ValueError("fs_hz must be positive")
    h = [_real(value, "coefficient") for value in coefficients]
    if not h:
        raise ValueError("Coefficients must be nonempty")
    result = []
    for frequency in frequencies_hz:
        f = _real(frequency, "frequency")
        if not 0 <= f <= fs / 2:
            raise ValueError("Frequency must be between DC and Nyquist")
        result.append(sum(coefficient * cmath.exp(-2j * math.pi * f * n / fs)
                          for n, coefficient in enumerate(h)))
    return result
