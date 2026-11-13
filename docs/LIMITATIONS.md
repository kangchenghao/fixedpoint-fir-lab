# Limits and interpretation

- Two constructed synthetic records are examples, not a representative EEG dataset or a training/holdout study. Microvolt labels specify a convenient numerical scale only.
- The small record is explicitly rounded to a quarter-uV grid by its generator. Its zero input-quantization error is deliberate and must not be generalized.
- The float reference is the same finite Hamming FIR design before quantization. It is not an ideal bandpass, continuous ground truth or a clinical gold standard.
- The 1% NRMSE line is illustrative. There is no clinical or downstream application acceptance claim.
- `fir_float` is intended for the documented finite signal range. Extremely large finite floating operands can overflow even though each input is finite.
- Python integers are arbitrary precision. The signed accumulator guard rejects overflow rather than modeling wrap or saturation. No target compiler, RTL, DSP, FPGA, area, clock, power or wall-time speedup is evaluated.
- The conservative bound covers direct terms and paired terms. The paired input sum can exceed the original 14-bit input range; a hardware implementation would need to size it separately.
- Direct/symmetric equality is algebraic and checked numerically; this is not a claim that one Python loop is faster. All cases include startup and tail.
- Input and coefficient rounding use binary-floating Python `round`; exact halfway representable test cases exercise ties to even. Near-boundary decimal inputs can have binary representation effects.
- `verify_integer_output` checks the mathematical output with a different loop organization. It does not prove accumulator width, optimizer correctness, all implementations or all future configurations.
- Tests and a written range argument are not Coq/Lean proofs. The suite includes deliberately broken cases to demonstrate its checking scope.
- Reproducing random numbers alone does not guarantee bit-identical transcendental functions on every platform. The report records per-run hashes, Python version and platform; numerical tolerances are used where appropriate.

Before extending this work, specify the new domain, rounding, overflow and observation contract. Recompute bounds and retain adverse cases, rather than selecting only examples that favor one precision.
