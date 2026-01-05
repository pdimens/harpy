def trunc_digits(x: float,y: int) -> float:
  '''Trucate the input float `x` at decimal digit `y` without rounding'''
  return float(f"%.{y}f" % x)

