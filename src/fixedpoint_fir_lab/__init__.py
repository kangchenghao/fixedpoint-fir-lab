"""Numerical FIR experiments; no EEG data or hardware model."""

from .core import (
    accumulator_bound, design_bandpass, fir_accumulators, fir_float, fir_integer,
    frequency_response, guard_accumulator, quantize_coefficients, quantize_input,
    verify_integer_output,
)

__version__ = "0.1.0"
__all__ = [
    "accumulator_bound", "design_bandpass", "fir_accumulators", "fir_float",
    "fir_integer", "frequency_response", "guard_accumulator",
    "quantize_coefficients", "quantize_input", "verify_integer_output",
]
