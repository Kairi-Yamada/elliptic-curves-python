"""Elliptic_Curve_over_Fp.py の動作を確認するテスト。

元のプログラムとこのファイルを同じフォルダーに置いて実行する。
    python -m unittest -v test_elliptic_curve

このファイルを直接実行することもできる。
    python test_elliptic_curve.py -v

Python 3.8 以上、標準ライブラリのみを使用する。
具体的な計算結果・例外処理・小さい有限群の性質を確認する。
"""

import unittest

from Elliptic_Curve_over_Fp import EllipticCurve, O, is_prime


class TestIsPrime(unittest.TestCase):
    """素数と、それ以外の整数を正しく区別する。"""

    def test_primes(self):
        for n in (2, 3, 5, 7, 11, 17, 97):
            with self.subTest(n=n):
                self.assertTrue(is_prime(n))

    def test_non_primes(self):
        # 平方数も含め、試し割りの上限付近の誤りを検出する。
        for n in (-7, -1, 0, 1, 4, 9, 25, 49, 91, 121):
            with self.subTest(n=n):
                self.assertFalse(is_prime(n))


class TestCurveCreation(unittest.TestCase):
    """係数の正規化と、楕円曲線の定義に必要な入力検査。"""

    def test_coefficients_are_normalized(self):
        curve = EllipticCurve(-15, 36, 17)
        self.assertEqual(curve.A, 2)
        self.assertEqual(curve.B, 2)
        self.assertEqual(curve.p, 17)

    def test_characteristic_below_five_is_rejected(self):
        for p in (-7, 0, 1, 2, 3, 4):
            with self.subTest(p=p):
                with self.assertRaises(ValueError):
                    EllipticCurve(2, 2, p)

    def test_composite_modulus_is_rejected(self):
        for p in (6, 9, 25, 49):
            with self.subTest(p=p):
                with self.assertRaises(ValueError):
                    EllipticCurve(2, 2, p)

    def test_singular_curves_are_rejected(self):
        # (A, B, p) = (2, 2, 5) では 4A^3 + 27B^2 = 140 ≡ 0。
        for a, b, p in ((0, 0, 5), (2, 2, 5), (0, 0, 97)):
            with self.subTest(A=a, B=b, p=p):
                with self.assertRaises(ValueError):
                    EllipticCurve(a, b, p)

    def test_non_integer_parameters_are_rejected(self):
        for position, name in enumerate(("A", "B", "p")):
            for value in (True, False, 1.5, "2", None):
                with self.subTest(parameter=name, value=value):
                    parameters = [2, 2, 17]
                    parameters[position] = value
                    with self.assertRaises(TypeError):
                        EllipticCurve(*parameters)


class TestPoints(unittest.TestCase):
    """点の生成、曲線上の判定、座標の正規化。"""

    def setUp(self):
        self.curve = EllipticCurve(2, 2, 17)
        self.P = (5, 1)

    def test_create_point(self):
        self.assertEqual(self.curve.point(5, 1), self.P)

    def test_coordinates_are_normalized(self):
        self.assertEqual(self.curve.point(22, -16), self.P)
        self.assertEqual(self.curve.normalize_point((22, -16)), self.P)

    def test_points_on_curve_are_recognized(self):
        for point in (self.P, (22, -16), O):
            with self.subTest(point=point):
                self.assertTrue(self.curve.is_on_curve(point))

    def test_point_outside_curve_is_recognized(self):
        # (0, 0) は y^2 = x^3 + 2x + 2 (mod 17) を満たさない。
        self.assertFalse(self.curve.is_on_curve((0, 0)))

    def test_point_outside_curve_is_rejected(self):
        with self.assertRaises(ValueError):
            self.curve.point(0, 0)
        with self.assertRaises(ValueError):
            self.curve.normalize_point((0, 0))

    def test_non_integer_coordinates_are_rejected(self):
        for value in (True, False, 1.5, "5", None):
            for point in ((value, 1), (5, value)):
                with self.subTest(point=point):
                    with self.assertRaises(TypeError):
                        self.curve.point(*point)
                    with self.assertRaises(TypeError):
                        self.curve.is_on_curve(point)

    def test_infinity_is_preserved_by_normalization(self):
        self.assertIs(self.curve.normalize_point(O), O)


class TestAdditionAndNegation(unittest.TestCase):
    """通常の加算、2倍、無限遠点、逆元を確認する。"""

    def setUp(self):
        self.curve = EllipticCurve(2, 2, 17)
        self.P = (5, 1)
        self.Q = (6, 3)

    def test_infinity_is_identity(self):
        self.assertEqual(self.curve.add(self.P, O), self.P)
        self.assertEqual(self.curve.add(O, self.P), self.P)
        self.assertIs(self.curve.add(O, O), O)

    def test_add_distinct_points(self):
        # 傾きは (3 - 1)/(6 - 5) = 2。
        # x3 = 2^2 - 5 - 6 ≡ 10, y3 = 2(5 - 10) - 1 ≡ 6。
        self.assertEqual(self.curve.add(self.P, self.Q), (10, 6))

    def test_double_point(self):
        # 2^(-1) ≡ 9, 傾きは (3・5^2 + 2)・9 ≡ 13 (mod 17)。
        self.assertEqual(self.curve.add(self.P, self.P), self.Q)

    def test_add_inverse_points(self):
        self.assertIs(self.curve.add(self.P, (5, 16)), O)

    def test_double_point_with_zero_y(self):
        curve = EllipticCurve(1, 0, 5)
        point = (0, 0)
        self.assertTrue(curve.is_on_curve(point))
        self.assertIs(curve.add(point, point), O)

    def test_addition_normalizes_input_before_comparison(self):
        # 整数の組としては異なるが、(22, -16) と (5, 1) は同じ点。
        self.assertEqual(self.curve.add((22, -16), self.P), self.Q)

    def test_negation(self):
        self.assertEqual(self.curve.neg(self.P), (5, 16))
        self.assertEqual(self.curve.neg((22, -16)), (5, 16))
        self.assertEqual(self.curve.neg((5, 16)), self.P)

    def test_negation_of_infinity(self):
        self.assertIs(self.curve.neg(O), O)

    def test_arithmetic_rejects_points_outside_curve(self):
        for left, right in (
            ((0, 0), self.P),
            (self.P, (0, 0)),
            ((0, 0), O),
            (O, (0, 0)),
        ):
            with self.subTest(P=left, Q=right):
                with self.assertRaises(ValueError):
                    self.curve.add(left, right)
        with self.assertRaises(ValueError):
            self.curve.neg((0, 0))


class TestScalarMultiplication(unittest.TestCase):
    """既知の整数倍と、負の整数・0・無限遠点を確認する。"""

    def setUp(self):
        self.curve = EllipticCurve(2, 2, 17)
        self.P = (5, 1)

    def test_known_positive_multiples(self):
        # E: y^2 = x^3 + 2x + 2 (mod 17), P = (5, 1)。
        # 期待値は固定した表を使い、mul や add の結果から作らない。
        expected = {
            1: (5, 1),
            2: (6, 3),
            3: (10, 6),
            4: (3, 1),
            5: (9, 16),
            6: (16, 13),
            7: (0, 6),
            8: (13, 7),
            9: (7, 6),
            10: (7, 11),
            11: (13, 10),
            12: (0, 11),
            13: (16, 4),
            14: (9, 1),
            15: (3, 16),
            16: (10, 11),
            17: (6, 14),
            18: (5, 16),
            19: O,
            20: (5, 1),
            37: (5, 16),
            38: O,
            1000: (0, 11),
        }
        for n, point in expected.items():
            with self.subTest(n=n):
                self.assertEqual(self.curve.mul(n, self.P), point)

    def test_negative_multiples(self):
        expected = {
            -1: (5, 16),
            -2: (6, 14),
            -3: (10, 11),
            -19: O,
            -20: (5, 16),
            -1000: (0, 6),
        }
        for n, point in expected.items():
            with self.subTest(n=n):
                self.assertEqual(self.curve.mul(n, self.P), point)

    def test_zero_multiple(self):
        self.assertIs(self.curve.mul(0, self.P), O)

    def test_multiples_of_infinity(self):
        for n in (-17, 0, 1, 38, 1000):
            with self.subTest(n=n):
                self.assertIs(self.curve.mul(n, O), O)

    def test_scalar_multiplication_normalizes_input(self):
        self.assertEqual(self.curve.mul(3, (22, -16)), (10, 6))

    def test_non_integer_scalars_are_rejected(self):
        for n in (True, False, 1.5, "2", None):
            with self.subTest(n=n):
                with self.assertRaises(TypeError):
                    self.curve.mul(n, self.P)

    def test_point_outside_curve_is_rejected(self):
        for n in (0, 1, -1):
            with self.subTest(n=n):
                with self.assertRaises(ValueError):
                    self.curve.mul(n, (0, 0))


class TestGroupLaws(unittest.TestCase):
    """3本の小さい曲線で、点を全列挙して群の性質を確認する。"""

    def setUp(self):
        self.examples = []
        for a, b, p in ((1, 0, 5), (0, 1, 7), (1, 6, 11)):
            curve = EllipticCurve(a, b, p)
            # is_on_curve 自体に誤りがあっても検出できるよう、
            # 点の列挙には曲線の方程式を直接使う。
            points = [O] + [
                (x, y)
                for x in range(p)
                for y in range(p)
                if (y * y - x * x * x - a * x - b) % p == 0
            ]
            self.examples.append((curve, points))

    def test_closure_and_commutativity(self):
        for curve, points in self.examples:
            point_set = set(points)
            for P in points:
                for Q in points:
                    with self.subTest(
                        A=curve.A, B=curve.B, p=curve.p, P=P, Q=Q
                    ):
                        result = curve.add(P, Q)
                        self.assertIn(result, point_set)
                        self.assertEqual(result, curve.add(Q, P))

    def test_associativity(self):
        for curve, points in self.examples:
            for P in points:
                for Q in points:
                    for R in points:
                        with self.subTest(
                            A=curve.A, B=curve.B, p=curve.p,
                            P=P, Q=Q, R=R
                        ):
                            left = curve.add(curve.add(P, Q), R)
                            right = curve.add(P, curve.add(Q, R))
                            self.assertEqual(left, right)


if __name__ == "__main__":
    unittest.main()
