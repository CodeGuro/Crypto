# library to aid in the math of twisted edwards curves
import random
from _sha3 import sha3_512


def int_to_bytes(n):
    """
    :param n: the integer to convert to bytes
    :return: bytes that represents the int in little endian
    """
    return int.to_bytes(n, int.bit_length(n), 'little')


def int_from_bytes(b):
    """
    :param b: the bytes in little endian to convert to an integer
    :return: the integer
    """
    return int.from_bytes(b, 'little')


def point_to_bytes(xy):
    """
    :param xy: the x,y tuple that represents the point
    :return: bytes in that represents the point
    """
    x, y = xy
    return int_to_bytes(x) + int_to_bytes(y)


def modInverse(a, p):
    """ Use Fermatt's little theorem to calculate module inverse modulus some prime p
    :param a: the constant to invert
    :param p: the modulus constant
    :return: b such that a^(p-2) = b (mod p)
    """
    return pow(a, p-2, p)


def legendre(a, p):
    return pow(a, (p - 1) // 2, p)


def mod_sqrt(n, p):
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


def twistEdPointAdd(p1, p2, a, d, p):
    """Point addition for twisted edwards curves
    :param p1: x1,y1 tuple for the first point
    :param p2: x2,y2 tuple for the second point
    :param a: the constant 'a' as used for the twisted edwards curve
    :param d: the constant 'd' as used for the twisted edwards curve
    :param p: the prime field p as used for the prime field
    :return: x3,y3 for the resulting point
    """

    if p1 == (None, None) and p2 == (None, None):
        raise Exception('adding 2 infi points')

    if p1 == (None, None):
        return p2

    if p2 == (None, None):
        return p1

    x1, y1 = p1
    x2, y2 = p2

    x3 = ((x1*y2 + x2*y1) % p * modInverse(1 + d*x1*x2*y1*y2, p)) % p
    y3 = ((y1*y2 - a*x1*x2) % p * modInverse(1 - d*x1*x2*y1*y2, p)) % p

    return x3, y3


def twistEdPointDouble(point, a, d, p):
    """Point doubling for twisted edwards curves
    :param point: x1, y1 tuple for the incoming point
    :param a: the constant 'a' as used for the twisted edwards curve
    :param d: the constant 'd' as used for the twisted edwards curve
    :param p: the prime field p as used for the prime field
    :return: x3, y3 for the resulting point
    """
    x1, y1 = point

    x3 = ((2*y1*x1 % p) * modInverse((a*x1*x1 % p) + (y1*y1 % p), p)) % p
    y3 = ((((y1*y1 % p) - (a*x1*x1 % p)) % p) * modInverse(2 - (a*x1*x1 % p) - (y1*y1 % p), p)) % p

    return x3, y3


def twistEdPointMultiply(point, scalar, a, d, p):
    """Point Multiplications for twisted edwards curves
    :param point: x1, y1 tuple for the incoming point
    :param scalar: The scalar n to multiply the point such that Q = nP
    :param a: the constant 'a' as used for the twisted edwards curve
    :param d: the constant 'd' as used for the twisted edwards curve
    :param p: the prime field p as used for the prime field
    :return: x3, y3 for the resulting point that is Q = nP
    """

    if scalar == 0:
        return None, None

    x1, y1 = point

    str_scalar_bin = str(bin(scalar))[2:]
    qx, qy = x1, y1

    for i in str_scalar_bin[1:]:
        qx, qy = twistEdPointDouble((qx, qy), a, d, p)
        if i == '1':
            qx, qy = twistEdPointAdd((qx, qy), (x1, y1), a, d, p)

    return qx, qy


def twistEdSign(gen_point, priv, hash, a, d, p, n):
    """Generate signatures for twisted edwards curves
    :param gen_point: x1, y1 of the generator point
    :param priv: the private key to sign with
    :param hash: hash of the message to sign
    :param a: the constant 'a' as used for the twisted edwards curve
    :param d: the constant 'd' as used for the twisted edwards curve
    :param p: the prime field p as used for the prime field
    :param n: the order of the generator point (order of the subgroup generated by gen_point)
    :return: (r, s) that is the signature
    """

    r = 0
    s = 0

    while r == 0 or s == 0:
        k = random.randrange(1, n-1)
        x1_q, y1_q = twistEdPointMultiply(gen_point, k, a, d, p)

        # r can be set arbitrarily, but we choose r to be x_q such that Q = uH + vG where H is the public key
        # This way we don't have to use more memory to pass more data for point Q which the verifier also needs.
        # The way the signature works: Q = uH + vG = kG where G is the the generator and k is the (temporary) priv key
        # k = inverse(s) * (z + r*d_a)  where z is the message and d_a is the (permanent) priv key
        # Q = uH + vG
        # Q = inverse(s)*rH +  inverse(s)*zG
        # Q = inverse(s)*r*d_a*G + inverse(s)*zG
        # Q = inverse(s)*(r*d_a + z)G
        # Q = kG
        # This allows us to give distinct values of r,s without giving any information about k or d_a private keys
        # as long as we don't re-use 'k'.
        r = x1_q
        s = modInverse(k, n) * (hash + r * priv) % n

    return r, s


def twistEdVerify(gen_point, pub_key, sig, hash, a, d, p, n):
    """Verify a signature
    :param gen_point: x, y of the generator point
    :param pub_key:
    :param sig: (r, s) that is the signature
    :param hash: the hash of the message to sign
    :param a: the constant 'a' as used for the twisted edwards curve
    :param d: the constant 'd' as used for the twisted edwards curve
    :param p: the prime field p as used for the prime field
    :param n: the order of the generator point (order of the subgroup generated by gen_point)
    :return: True if the signature is valid, False if not.
    """
    r, s = sig

    s_inv = modInverse(s, n)
    u1 = s_inv * r % n
    u2 = s_inv * hash % n
    px, py = twistEdPointAdd(twistEdPointMultiply(pub_key, u1, a, d, p),
                             twistEdPointMultiply(gen_point, u2, a, d, p),
                             a, d, p)

    return r % n == px % n


def twistEdSign2(gen_point, priv, hash, a, d, p, n):
    """A safer implementation of twistEdSign, using a different scheme but the same principles as the previous implementation
    :param gen_point: x1, y1 of the generator point
    :param priv: the private key to sign with
    :param hash: hash of the message to sign
    :param a: the constant 'a' as used for the twisted edwards curve
    :param d: the constant 'd' as used for the twisted edwards curve
    :param p: the prime field p as used for the prime field
    :param n: the order of the generator point (order of the subgroup generated by gen_point)
    :return: (R, s) that is the signature
    """

    # The idea this time is we use S = R + hP, where 'S' and 'R' are given via the signature and P is the public key.
    # 'h' is the hash(P + R + z) where 'z' is the message (or message hash)
    # The client has full disclosure on constructing 'h' due to the given parameters {P, R, z}
    # This yields the equation sG = rG + hpG, where rG = R and pG = P
    # Note: this means 'p' must be kept entirely secret, which it is due to the elliptic curve discrete logarithm problem
    # To keep the signature rigid, we make h = hash(hash(p)+bytes(R)+bytes(P)). This is a hard problem for someone trying to make a bogus signature
    # because the malicious entity would need to know private key (p) in order to calculate a bogus h', r'
    # via modular inverse and use the equation v' = h*v*inv(v)*bogus_v to a malicious signature: s'G = u'G + hv'G = Q' = R' + hP'
    # That (P') is essentially the ECDLP due to hP = h(vG). Ha Ha! You~ Can't~ Get~ it~! Na na na na na~!


    P = twistEdPointMultiply(gen_point, priv, a, d, p)
    bytes_msg = int_to_bytes(hash)
    r = int_from_bytes(sha3_512(sha3_512(int_to_bytes(priv)).digest() + bytes_msg).digest()) % n
    R = twistEdPointMultiply(gen_point, r, a, d, p)
    bytes_R = sha3_512(point_to_bytes(R)).digest()
    bytes_P = sha3_512(point_to_bytes(P)).digest()
    h = int_from_bytes(sha3_512(bytes_R + bytes_P + bytes_msg).digest()) % n
    s = (r + h*priv % n) % n

    return r, s  # alternatively, return R, s


def twistEdVerify2(gen_point, pub_key, sig, hash, a, d, p, n):
    """Verify a signature, using an implementation compatible with twistEdSign2
    :param gen_point: x, y of the generator point
    :param pub_key:
    :param sig: (r, s) (or (u,s)) that is the signature
    :param hash: the hash of the message to sign
    :param a: the constant 'a' as used for the twisted edwards curve
    :param d: the constant 'd' as used for the twisted edwards curve
    :param p: the prime field p as used for the prime field
    :param n: the order of the generator point (order of the subgroup generated by gen_point)
    :return: True if the signature is valid, False if not.
    """
    u, s = sig
    P = pub_key
    sG = twistEdPointMultiply(gen_point, s, a, d, p)
    R = twistEdPointMultiply(gen_point, u, a, d, p)
    bytes_R = sha3_512(point_to_bytes(R)).digest()
    bytes_P = sha3_512(point_to_bytes(P)).digest()
    bytes_msg = int_to_bytes(hash)
    h = int_from_bytes(sha3_512(bytes_R + bytes_P + bytes_msg).digest()) % n
    hP = twistEdPointMultiply(pub_key, h, a, d, p)

    _sG = twistEdPointAdd(R, hP, a, d, p)

    valid = sG == _sG
    return valid


def ECDH(rem_pub_key, priv_key, a, d, p):
    """
    :param rem_pub_key: the remote public key (the client should not have the private key for this)
    :param priv_key: the local private key (the remote should not have this private key)
    :param a: the constant 'a' as used for the twisted edwards curve
    :param d: the constant 'd' as used for the twisted edwards curve
    :param p: the prime field p as used for the prime field
    :return:
    """
    shared_secret = twistEdPointMultiply(rem_pub_key, priv_key, a, d, p)
    return point_to_bytes(shared_secret)




# E-521 curves for 519 bits of security because screw NIST! Tell us where the unexplained seed was for the hash dammit!
message = b'this is some message'
hash_sample = int_from_bytes(message)
generator_point = (0x752cb45c48648b189df90cb2296b2878a3bfd9f42fc6c818ec8bf3c9c0c6203913f6ecc5ccc72434b1ae949d568fc99c6059d0fb13364838aa302a940a2f19ba6c, 0xc)
gen_order = 0x7ffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffd15b6c64746fc85f736b8af5e7ec53f04fbd8c4569a8f1f4540ea2435f5180d6b
prime_field = 0x1ffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff
edD = -376014
edA = 1

print('Generating key pair for Alice')
priv_key = int_from_bytes(b'alicealicealicealicealicealiceal') % gen_order
pub_key = twistEdPointMultiply(generator_point, priv_key, edA, edD, prime_field)

print('Sigining message ver1')
sig = twistEdSign(generator_point, priv_key, hash_sample, edA, edD, prime_field, gen_order)

print('Verifying message ver1')
ver = twistEdVerify(generator_point, pub_key, sig, hash_sample, edA, edD, prime_field, gen_order)

print('Signing message ver2')
sig2 = twistEdSign2(generator_point, priv_key, hash_sample, edA, edD, prime_field, gen_order)

print('Verifying message ver2')
ver2 = twistEdVerify2(generator_point, pub_key, sig2, hash_sample, edA, edD, prime_field, gen_order)

print('validity ver1 : ver2')
print(str(ver) + ' : ' + str(ver2))

print('Generating key pair for Bob')
priv_key2 = int_from_bytes(b'bobbobbobbobbobbo') % gen_order
pub_key2 = twistEdPointMultiply(generator_point, priv_key2, edA, edD, prime_field)

print('Performing ECHD handshake')
sec1 = ECDH(pub_key, priv_key2, edA, edD, prime_field)
sec2 = ECDH(pub_key2, priv_key, edA, edD, prime_field)

print('ECDH verification: ' + str(sec1 == sec2))
