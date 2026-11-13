import math
import random
import unittest

from fixedpoint_fir_lab import (
    accumulator_bound, design_bandpass, fir_accumulators, fir_float, fir_integer,
    frequency_response, guard_accumulator, quantize_coefficients, quantize_input,
    verify_integer_output,
)


class DesignTests(unittest.TestCase):
    def test_symmetry_and_center_normalization(self):
        taps = design_bandpass()
        self.assertEqual(len(taps), 31)
        self.assertEqual(taps, taps[::-1])
        self.assertAlmostEqual(abs(frequency_response(taps, [19])[0]), 1.0, places=13)

    def test_transition_centers_and_stopband(self):
        dc, low, high, stop = map(abs, frequency_response(design_bandpass(), [0, 8, 30, 55]))
        self.assertLess(dc, 0.03)
        self.assertTrue(0.48 < low < 0.52)
        self.assertTrue(0.48 < high < 0.52)
        self.assertLess(stop, 0.001)

    def test_linear_phase_where_response_is_nonzero(self):
        for frequency, response in zip((8, 12, 19, 30), frequency_response(design_bandpass(), [8, 12, 19, 30])):
            # Remove the independently known 15-sample linear phase term.
            compensated = response * complex(math.cos(2 * math.pi * frequency * 15 / 160),
                                               math.sin(2 * math.pi * frequency * 15 / 160))
            self.assertAlmostEqual(compensated.imag, 0.0, places=12)

    def test_invalid_design_and_response(self):
        for args in ((0, 8, 30, 31), (160, 0, 30, 31), (160, 30, 8, 31),
                     (160, 8, 80, 31), (160, 8, 30, 30), (160, 8, 30, True),
                     (float("nan"), 8, 30, 31)):
            with self.subTest(args=args), self.assertRaises(ValueError):
                design_bandpass(*args)
        for frequency in (-1, 81, float("nan")):
            with self.subTest(frequency=frequency), self.assertRaises(ValueError):
                frequency_response([1], [frequency])


class QuantizationTests(unittest.TestCase):
    def test_input_ties_to_even_both_signs(self):
        self.assertEqual(quantize_input([0.125, 0.375, -0.125, -0.375, 0.625, -0.625]), [0, 2, 0, -2, 2, -2])

    def test_coefficient_ties_to_even_both_signs(self):
        self.assertEqual(quantize_coefficients([0.5 / 32768, 1.5 / 32768, -0.5 / 32768, -1.5 / 32768]), [0, 2, 0, -2])

    def test_input_domain_and_no_silent_clipping(self):
        self.assertEqual(quantize_input([-2048, 2047.75]), [-8192, 8191])
        for value in (-2048.25, 2048, float("nan"), float("inf"), True):
            with self.subTest(value=value), self.assertRaises(ValueError):
                quantize_input([value])

    def test_coefficient_domain(self):
        self.assertEqual(quantize_coefficients([-1, 32767 / 32768]), [-32768, 32767])
        for value in (1, -1.1, float("nan")):
            with self.subTest(value=value), self.assertRaises(ValueError):
                quantize_coefficients([value])

    def test_negative_output_uses_floor_not_truncation(self):
        actual = fir_integer([-1, 1], [1], output_right_shift=1)
        self.assertEqual(actual, [-1, 0])
        with self.assertRaises(ValueError):
            verify_integer_output([-1, 1], [1], [0, 0], output_right_shift=1)

    def test_output_scales_represent_the_same_uv(self):
        for shift, units in ((15, 4), (13, 16), (11, 64)):
            result = fir_integer([4], [16384], output_right_shift=shift)
            self.assertEqual(result[0] / units, 0.5)


class ConvolutionTests(unittest.TestCase):
    def test_hand_computed_negative_warmup_and_tail(self):
        self.assertEqual(fir_accumulators([3, -2], [1, 2, 1]), [3, 4, -1, -2])
        self.assertEqual(fir_float([3, -2], [1, 2, 1]), [3.0, 4.0, -1.0, -2.0])

    def test_impulse_is_the_coefficient_sequence(self):
        taps = design_bandpass()
        q15 = quantize_coefficients(taps)
        self.assertEqual(fir_accumulators([1], q15), q15)
        self.assertEqual(fir_float([1], taps), taps)

    def test_prefix_is_causal_even_with_tail(self):
        taps = quantize_coefficients(design_bandpass())
        prefix = fir_accumulators([10, -17, 25], taps)
        extended = fir_accumulators([10, -17, 25, -8192, 8191], taps)
        self.assertEqual(prefix[:3], extended[:3])
        self.assertEqual(len(extended), 5 + 31 - 1)

    def test_symmetric_extremes_and_seeded_inputs(self):
        rng = random.Random(410)
        q15 = quantize_coefficients(design_bandpass())
        for x in ([0], [-8192] * 40, [8191] * 40, [-8192, 8191] * 20,
                  [rng.randrange(-8192, 8192) for _ in range(71)]):
            with self.subTest(length=len(x), first=x[0]):
                direct = fir_accumulators(x, q15, accumulator_bits=32)
                self.assertEqual(direct, fir_accumulators(x, q15, method="symmetric", accumulator_bits=32))
                for shift in (11, 13, 15):
                    output = [value >> shift for value in direct]
                    self.assertEqual(verify_integer_output(x, q15, output, output_right_shift=shift), len(output))

    def test_tampered_value_and_missing_tail_are_rejected(self):
        x, h = [3, -2], [1, 2, 1]
        expected = [3, 4, -1, -2]
        for wrong in ([3, 4, 0, -2], expected[:-1], [3, 4, -1, -2, 0]):
            with self.subTest(output=wrong), self.assertRaises(ValueError):
                verify_integer_output(x, h, wrong, output_right_shift=0)

    def test_invalid_domains_and_methods(self):
        for x, h, kwargs in (([], [1], {}), ([1], [], {}), ([8192], [1], {}),
                             ([True], [1], {}), ([1], [32768], {}),
                             ([1], [1], {"method": "unknown"}),
                             ([1], [1, 2], {"method": "symmetric"}),
                             ([1], [1, 2, 3], {"method": "symmetric"}),
                             ([1], [1], {"output_right_shift": -1}),
                             ([1], [1], {"output_right_shift": 32})):
            with self.subTest(x=x, h=h, kwargs=kwargs), self.assertRaises(ValueError):
                fir_integer(x, h, **kwargs)


class AccumulatorTests(unittest.TestCase):
    def test_bound_covers_all_partial_sums(self):
        q15 = quantize_coefficients(design_bandpass())
        bound = accumulator_bound(q15)
        self.assertEqual(bound, 389455872)
        self.assertLessEqual(bound, 2**31 - 1)
        self.assertEqual(guard_accumulator(-(2**31)), -(2**31))
        self.assertEqual(guard_accumulator(2**31 - 1), 2**31 - 1)

    def test_positive_and_negative_overflow_are_rejected(self):
        for value in (2**31, -(2**31) - 1):
            with self.subTest(value=value), self.assertRaises(OverflowError):
                guard_accumulator(value)
        for sample in (8191, -8192):
            for method in ("direct", "symmetric"):
                with self.subTest(sample=sample, method=method), self.assertRaises(OverflowError):
                    fir_accumulators([sample] * 33, [32767] * 33, method=method, accumulator_bits=32)
        # Python reference remains defined: overflow is not silently wrapped.
        self.assertGreater(max(fir_accumulators([8191] * 33, [32767] * 33)), 2**31 - 1)


if __name__ == "__main__":
    unittest.main()
