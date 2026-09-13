import sys

PARAMS = [(2, 2), (2, 3), (3, 2), (3, 3), (4, 2)]
MODULUS_BOUND = 64


def is_prime(n):
    if n < 2:
        return False
    d = 2
    while d * d <= n:
        if n % d == 0:
            return False
        d += 1
    return True


def cios(a, b, m, w, T):
    mask = (1 << w) - 1
    mp = (-pow(m, -1, 1 << w)) % (1 << w)
    A = [0] * (T + 2)
    aw = [(a >> (w * i)) & mask for i in range(T)]
    bw = [(b >> (w * i)) & mask for i in range(T)]
    mw = [(m >> (w * i)) & mask for i in range(T)]
    for i in range(T):
        C = 0
        for j in range(T):
            t = A[j] + aw[j] * bw[i] + C
            A[j] = t & mask
            C = t >> w
        t = A[T] + C
        A[T] = t & mask
        A[T + 1] = t >> w
        mm = (A[0] * mp) & mask
        C = (A[0] + mm * mw[0]) >> w
        for j in range(1, T):
            t = A[j] + mm * mw[j] + C
            A[j - 1] = t & mask
            C = t >> w
        t = A[T] + C
        A[T - 1] = t & mask
        A[T] = A[T + 1] + (t >> w)
    return sum(A[i] << (w * i) for i in range(T + 1))


def main():
    evaluations = 0
    hits_prime = 0
    hits_composite = 0
    for w, T in PARAMS:
        hi = min(1 << (w * T), MODULUS_BOUND)
        for m in range(3, hi, 2):
            prime = is_prime(m)
            for a in range(m):
                for b in range(m):
                    evaluations += 1
                    if cios(a, b, m, w, T) == m:
                        if prime:
                            hits_prime += 1
                        else:
                            hits_composite += 1

    print("CIOS Montgomery post-condition reachability")
    print("  (w,T) pairs        %s" % (PARAMS,))
    print("  moduli             all odd m, 3 <= m < min(2^(wT), %d)" % MODULUS_BOUND)
    print("  operands           all (a,b), 0 <= a,b < m")
    print("  evaluations        %d" % evaluations)
    print("  A==m occurrences   %d" % (hits_prime + hits_composite))
    print("    prime modulus    %d" % hits_prime)
    print("    composite        %d" % hits_composite)

    if hits_prime != 0:
        print("\nFAIL: A==m reached under a prime modulus")
        return 1
    if hits_composite == 0:
        print("\nFAIL: A==m never reached; sweep proves nothing")
        return 1
    print("\nok: A==m reachable only for composite moduli")
    return 0


if __name__ == "__main__":
    sys.exit(main())
