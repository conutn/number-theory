import math
import random
from math import e, pi

def primitive_roots(q):
    candidates = range(1, q)
    lst = []
    for n in candidates:
        s = {pow(n, p, q) for p in range(q - 1)}
        if len(s) == q - 1:
            lst.append(n)
    return lst

def v(n, g, q):
    if n % q == 0: return 0
    for i in range(0, q - 1):
        if pow(g, i, q) == n % q:
            return i

def dirichlet_character(i, g, q):
    omega = e ** (i * 2j * math.pi / (q - 1))
    return (lambda x: omega ** v(x, g, q) if x % q != 0 else 0)

def L(s, chi, lim):
    total = 0
    for i in range(1, lim + 1):
        total += chi(i) / i ** s
    return total

def legendre(a, q):
    if a % q == 0: return 0
    p = pow(a, (q-1)//2, q)
    if p == q - 1: return -1
    return 1

def G(q):
    total = 0
    for i in range(1, q):
        total += legendre(i, q) * e ** (i * 2j * math.pi / q)
    return total
    

q = int(input("Input an odd prime:"))
assert all(q%i for i in range(2, int(q**(1/2)) + 1)) and q > 2, "q is not an odd prime"

roots = primitive_roots(q)
print(f"Primitive roots mod {q}: {roots}")

g = random.choice(roots)
assert g in roots, "g is not a primitive root"
print(f"Using primitive root {g} mod {q}")

chi = dirichlet_character((q-1)//2, g, q)
print(f"L_-1(1) estimate: {L(1, chi, 1000)}")

if q % 4 == 3:
    s = sum(i * legendre(i, q) for i in range(1, q))
    print(f"True value: {-s*pi / q**(3/2)}")
if q % 4 == 1:
    s = sum(math.log(2*math.sin(math.pi*i/q)) * legendre(i, q) for i in range(1, q))
    print(f"True value: {-s / q**(1/2)}")
print("Gauss sum G(1) mod q:", G(q), "≈", (1 + 1j ** (-q)) / (1 + 1j ** (-1)) * q ** 0.5)
