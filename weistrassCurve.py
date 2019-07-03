
def legendre(a, p):
    return pow(a, (p - 1) // 2, p)


def tonelli(n, p):
    assert legendre(n, p) == 1, "not a square (mod p)"
    q = p - 1
    s = 0
    while q % 2 == 0:
        q //= 2
        s += 1
    if s == 1:
        return pow(n, (p + 1) // 4, p)
    for z in range(2, p):
        if p - 1 == legendre(z, p):
            break
    c = pow(z, q, p)
    r = pow(n, (q + 1) // 2, p)
    t = pow(n, q, p)
    m = s
    t2 = 0
    while (t - 1) % p != 0:
        t2 = (t * t) % p
        for i in range(1, m):
            if (t2 - 1) % p == 0:
                break
            t2 = (t2 * t2) % p
        b = pow(c, 1 << (m - i - 1), p)
        r = (r * b) % p
        c = (b * b) % p
        t = (t * c) % p
        m = i
    return r


def modinv(a, n):  # Extended Euclidean Algorithm/'division' in elliptic curves
    lm, hm = 1, 0
    low, high = a % n, n
    while low > 1:
        ratio = high // low
        nm, new = hm - lm * ratio, high - low * ratio
        lm, low, hm, high = nm, new, lm, low
    return lm % n


def find_factors(k):
    """ Find the prime factors of an integer
    :param k: int
    :return: list of prime factors
    """
    factors = list()

    for i in range(1, k+1):
        if 0 == k % i:
            factors.append(i)

    return factors


def ECadd(xp, yp, xq, yq, Fp):  # Not true addition, invented for EC. It adds Point-P with Point-Q.
    if xp == xq and yp == -yq % Fp:
        return None, None
    if xp is None and yp is None and xq is not None and yq is not None:
        return xq, yq
    if xp is not None and yp is not None and xq is None and yq is None:
        return xp, yp
    if xp is None and yp is None and xq is None and yq is None:
        raise Exception('Adding 2 infinity points')
    m = ((yq - yp) * modinv(xq - xp, Fp)) % Fp
    xr = (m * m - xp - xq) % Fp
    yr = (m * (xp - xr) - yp) % Fp
    return xr, yr


def ECdouble(xp, yp, Fp):  # EC point doubling,  invented for EC. It doubles Point-P.
    if xp is None and yp is None or yp == 0:
        return None, None
    LamNumer = 3 * xp * xp + Acurve
    LamDenom = 2 * yp
    Lam = (LamNumer * modinv(LamDenom, Fp)) % Fp
    xr = (Lam * Lam - 2 * xp) % Fp
    yr = (Lam * (xp - xr) - yp) % Fp
    return xr, yr


def EccMultiply(xs, ys, Scalar, Fp):  # Double & add. EC Multiplication, Not true multiplication
    if Scalar == 0:
        return None, None
    ScalarBin = str(bin(Scalar))[2:]
    Qx, Qy = xs, ys
    for i in range(1, len(ScalarBin)):  # This is invented EC multiplication.
        Qx, Qy = ECdouble(Qx, Qy, Fp)  # print "DUB", Qx; print
        if ScalarBin[i] == "1":
            Qx, Qy = ECadd(Qx, Qy, xs, ys, Fp)  # print "ADD", Qx; print
    return Qx, Qy


# Super simple Elliptic Curve Presentation. No imported libraries, wrappers, nothing. # For educational purposes only
N_Field_Prime = 0xffffffff00000001000000000000000000000000ffffffffffffffffffffffff  # The proven prime
Acurve = -3
Bcurve = 41058363725152142129326129780047268409114441015993725554835256314039467401291  # This defines the curve. y^2 = x^3 + Acurve * x + Bcurve

# ensure that the curve is not singular
assert(4*(Acurve**3) + 27*(Bcurve**2) != 0)

# Compute N = #E(F_p) that is the order (cardinality) of the curve over the prime field 'Pcurve'
from pyschoof.reduced_computation_schoof import reduced_computation_schoof_algorithm as schoof
n_group_order = schoof(N_Field_Prime, Acurve, Bcurve)  # This must be a prime number

# We need to choose a generator and the order of the generator. Take some point P in the curve.
Px = 6
Py2 = (Px**3 + Acurve*Px + Bcurve) % N_Field_Prime
Py = tonelli(Py2, N_Field_Prime)  # use the quadratic residue properties to calculate modular square root

# To find the order of the subgroups of E(F_p) we must count the factors #E(F_p) by that point
list_factors = find_factors(n_group_order)

# find the smallest group order and l such that lh = n_group_order
r = n_group_order
l = 5

# we now have a group order 'r'
cofactor = n_group_order // r

generator = EccMultiply(Px, Py, cofactor, N_Field_Prime)


test_point = ()
l_group_cardinality = 1
for f in list_factors:
    test_point = EccMultiply(Px, Py, f, N_Field_Prime)
    if test_point == (None, None):
        l_group_cardinality = f
        break

cofactor = n_group_order // l_group_cardinality
print('The group cardinality is ' + str(l_group_cardinality) + ', and the cofactor is: ' + str(cofactor))


# we will now try to find a subgroup with cardinality N // cofactor

other_point = EccMultiply(Px, Py, cofactor, N_Field_Prime)
other_group_cardinality = 1
for f in list_factors:
    other_group_cardinality = f
    test_point = EccMultiply(other_point[0], other_point[1], f, N_Field_Prime)
    if test_point == (None, None):
        break

print('other group cardinality is: ' + str(other_group_cardinality))



GPoint = EccMultiply(Px, Py, cofactor, N_Field_Prime)
Gx = GPoint[0]
Gy = GPoint[1]

# Individual Transaction/Personal Information
privKey = N_Field_Prime - 1  # replace with any private key
RandNum = 22  # replace with a truly random number
HashOfThingToSign = 6  # the hash of your message/transaction


print()
print("******* Public Key Generation *********")
xPublicKey, yPublicKey = EccMultiply(Gx, Gy, privKey, N_Field_Prime)
print("the private key (in base 10 format):")
print(privKey)
print()
print("the uncompressed public key (starts with '04' & is not the public address):")
print("04", xPublicKey, yPublicKey)

print()
print("******* Signature Generation *********")
xRandSignPoint, yRandSignPoint = EccMultiply(Gx, Gy, RandNum, N_Field_Prime)
r = xRandSignPoint % n_group_order
print("r =", r)
s = ((HashOfThingToSign + r * privKey) * (modinv(RandNum, n_group_order))) % n_group_order
print("s =", s)

print()
print("******* Signature Verification *********>>")
w = modinv(s, n_group_order)
xu1, yu1 = EccMultiply(Gx, Gy, (HashOfThingToSign * w) % n_group_order, N_Field_Prime)
xu2, yu2 = EccMultiply(xPublicKey, yPublicKey, (r * w) % n_group_order, N_Field_Prime)
x, y = ECadd(xu1, yu1, xu2, yu2, N_Field_Prime)
print(r == x)
print()
