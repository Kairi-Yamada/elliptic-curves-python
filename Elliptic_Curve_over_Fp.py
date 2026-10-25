from math import isqrt


def is_prime(n):
  if n <= 1:
    return False

  if n == 2:
    return True

  if n % 2 == 0:
    return False

  for i in range(3, isqrt(n) + 1, 2):
    if n % i == 0:
      return False

  return True


O = None


class EllipticCurve:
  def __init__(self, A, B, p):
    if isinstance(A, bool) or not isinstance(A, int):
      raise TypeError("A must be an integer")

    if isinstance(B, bool) or not isinstance(B, int):
      raise TypeError("B must be an integer")

    if isinstance(p, bool) or not isinstance(p, int):
      raise TypeError("p must be an integer")

    self.p = p

    if self.p < 5:
      raise ValueError("p must be at least 5")

    if not is_prime(self.p):
      raise ValueError("p must be prime number")

    self.A = A % self.p
    self.B = B % self.p

    if (4 * pow(self.A, 3, self.p) + 27 * pow(self.B, 2, self.p)) % self.p == 0:
      raise ValueError("singular curve")


  def point(self, x, y):
    if isinstance(x, bool) or not isinstance(x, int):
      raise TypeError("x must be an integer")
    
    if isinstance(y, bool) or not isinstance(y, int):
      raise TypeError("y must be an integer")
    
    P = (x % self.p, y % self.p)

    if not self.is_on_curve(P):
      raise ValueError("point is not on the curve")

    return P


  def is_on_curve(self, P):
    if P is O:
      return True

    x, y = P

    if isinstance(x, bool) or not isinstance(x, int):
      raise TypeError("x must be an integer")

    if isinstance(y, bool) or not isinstance(y, int):
      raise TypeError("y must be an integer")

    return pow(y, 2, self.p) == (pow(x, 3, self.p) + self.A * x + self.B) % self.p


  def normalize_point(self, P):
    if P is O:
      return O

    x, y = P

    return self.point(x, y)


  def add(self, P, Q):
    P = self.normalize_point(P)
    Q = self.normalize_point(Q)

    if P is O:
      return Q

    if Q is O:
      return P

    x1, y1 = P
    x2, y2 = Q

    if x1 == x2 and (y1 + y2) % self.p == 0:
      return O

    if P != Q:
      m = ((y2 - y1) * pow(x2 - x1, -1, self.p)) % self.p

    else:
      m = ((3 * pow(x1, 2, self.p) + self.A) * pow(2 * y1, -1, self.p)) % self.p

    x3 = (pow(m, 2, self.p) - x1 - x2) % self.p
    y3 = (m * (x1 - x3) - y1) % self.p

    R = self.point(x3, y3)

    return R


  def neg(self, P):
    P = self.normalize_point(P)

    if P is O:
      return O

    x, y = P

    Q = self.point(x, -y)

    return Q


  def mul(self, n, P):
    if isinstance(n, bool) or not isinstance(n, int):
      raise TypeError("n must be an integer")

    P = self.normalize_point(P)

    if n == 0 or P is O:
      return O

    if n < 0:
      return self.mul(-n, self.neg(P))

    exp = n
    res = O
    base = P

    while exp > 0:
      if exp % 2 == 1:
        res = self.add(res, base)

      base = self.add(base, base)
      exp = exp // 2

    return res