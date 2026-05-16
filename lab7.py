import math

# Parametry wejściowe z tabeli
p1, t1, vu1, delta1 = 1.0, 298.0, 0.74, 59.0
p2, t2, vu2, delta2 = 5.0, 298.0, 0.74, 31.0
p3, t3, vu3, delta3 = 1.0, 375.0, 1.00, 65.0

# Wyznaczenie laminarnej prędkości spalania (SL)
sl1 = vu1 * math.sin(math.radians(delta1 / 2))
sl2 = vu2 * math.sin(math.radians(delta2 / 2))
sl3 = vu3 * math.sin(math.radians(delta3 / 2))

# Obliczenie współczynnika ciśnieniowego beta (zmiana ciśnienia, stała temperatura)
beta = math.log(sl2 / sl1) / math.log(p2 / p1)

# Obliczenie współczynnika temperaturowego alpha (zmiana temperatury, stałe ciśnienie)
alpha = math.log(sl3 / sl1) / math.log(t3 / t1)

# Wypisanie wyników w czytelnym formacie
print(f"SL1 (ref): {sl1 * 100:.2f} cm/s")
print(f"SL2 (ciśnienie): {sl2 * 100:.2f} cm/s")
print(f"SL3 (temperatura): {sl3 * 100:.2f} cm/s")
print(f"Współczynnik beta: {beta:.3f}")
print(f"Współczynnik alpha: {alpha:.3f}")