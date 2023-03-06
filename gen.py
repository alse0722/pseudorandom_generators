import argparse


class Generator():
    @staticmethod
    def linear_congruential(nc, iv):
        values = []
        n = nc
        a, c, m, x0 = iv
        for _ in range(n):
            xn = (a*x0 + c) % m
            values.append(xn)
            x0 = xn

        return values

    @staticmethod
    def additive(nc, iv):
        values = []
        n = nc
        s = iv[0]
        c, m = iv[-2], iv[-1]
        a_s, xs, = iv[1:s+1], iv[s+1:2*s+1]
        for i in range(n):
            xi = 0
            for j in range(s):
                xi += a_s[j] * xs[i+j]
            xi += c
            xi %= m
            xs.append(xi)
            values.append(xi)

        return values

    @staticmethod
    def lfsr(nc, iv):
        values = []
        n = nc
        p, s, m = iv[0], iv[1], iv[2]
        js = iv[3:m+3]
        y1p = iv[-1]
        x_init = [int(d) for d in bin(y1p)[2:]]
        x_init = [0] * (p - len(x_init)) + x_init
        x_all = []
        for _ in range(n):
            for _ in range(p):
                next_bit = 0
                for jt in js:
                    next_bit ^= x_init[jt-1]
                if len(x_all) == s:
                    values.append(int(''.join(str(x) for x in x_all), 2))
                    x_all.clear()
                if len(values) == n:
                    break
                x_init.append(next_bit)
                x_all.append(next_bit)
                x_init.pop(0)

        return values

    @staticmethod
    def dopoln(str_key, p):
        if len(str_key) < p:
            for i in range(p - len(str_key)):
                str_key = '0' + str_key
        return str_key

    @staticmethod
    def _lfsr(n, p, s, m, j_mas, y):
        a = []

        for i in range(p):
            a.append(0)
        for i in range(m):
            a[int(j_mas[i]) - 1] = 1

        potok = ''
        str_key = format(y, 'b')
        str_key = Generator.dopoln(str_key, p)
        potok += str_key
        x = 0
        indicator = 0
        indicator2 = 0

        for i in range(s * n):
            x = 0
            indicator2 = indicator
            for j in range(p):
                if int(a[j]) == 1:
                    x += int(potok[indicator2])
                indicator2 += 1

            x = x % 2
            potok += str(x)
            indicator += 1

        res = []

        amount = (len(potok) - p) // s
        help = ''
        k = p
        for i in range(amount):
            help = ''
            for j in range(s):
                help += potok[k]
                k += 1
            res.append(int(help, 2))

        return res

    @staticmethod
    def five_parameters(nc, iv):
        values = []
        n = nc
        p, w = iv[0], iv[1]
        y1p = iv[-1]
        q1, q2, q3 = iv[2], iv[3], iv[4]
        x_init = [int(d) for d in bin(y1p)[2:]]
        x_init = [0] * (p - len(x_init)) + x_init
        x_all = x_init[:]
        for _ in range(n):
            for _ in range(w):
                next_bit = x_all[0] ^ x_all[q1] ^ x_all[q2] ^ x_all[q3]
                x_all.append(next_bit)
                x_all.pop(0)
            values.append(int(''.join(str(x)
                          for x in x_all[-1:-w-1:-1][::-1]), 2))
        return values

    @staticmethod
    def rsa(nc, iv):
        values = []
        c = nc
        n, e, x0, w, l = iv
        x_all = []
        for _ in range(c):
            x0 = pow(x0, e, n)
            xi_b = [int(d) for d in bin(x0)[2:]]
            xi_b = [0] * (w - len(xi_b)) + xi_b
            xi_b = xi_b[::-1][:w][::-1]
            x_all += xi_b
            values.append(int(''.join(str(x) for x in x_all[:l]), 2))
            x_all = x_all[l:]
        return values

    @staticmethod
    def bbs(nc, iv):
        values = []
        c = nc
        n, x0, l = iv
        for _ in range(c):
            x_bin = []
            for _ in range(l):
                x0 = pow(x0, 2, n)
                x_bin.append(x0 % 2)
            values.append(int(''.join(str(x) for x in x_bin), 2))
        return values

    @staticmethod
    def rc4(nc, iv):
        values = []
        n = nc
        w, k = iv[0], iv[1:]
        s = list(range(256))
        j = 0
        for i in range(256):
            j = (j + s[i] + k[i]) % 256
            s[i], s[j] = s[j], s[i]

        i = 0
        j = 0

        x_all = []
        t = (w * n)//8 + 1
        for _ in range(t):
            i = (i + 1) % 256
            j = (j + s[i]) % 256
            s[i], s[j] = s[j], s[i]
            x = s[(s[i] + s[j]) % 256]

            xi_b = [int(d) for d in bin(x)[2:]]
            xi_b = [0] * (8 - len(xi_b)) + xi_b
            # print(xi_b)
            x_all += xi_b

        xi_b = []
        for i in range(0, len(x_all)):
            xi_b.append(x_all[i])
            if (i + 1) % w == 0:
                values.append(int(''.join(str(x) for x in xi_b), 2))
                xi_b = []

        return values

    @staticmethod
    def nfsr(nc, iv):
        values = []
        c = nc
        k, w = iv[0], iv[1]
        p, x = [iv[2]], [iv[3]]
        m = [iv[4]]
        j = [iv[5:5+m[0]]]
        p.append(iv[5+m[0]])
        x.append(iv[6+m[0]])
        m.append(iv[7+m[0]])
        j.append(iv[8+m[0]:8+m[0]+m[1]])
        p.append(iv[8+m[0]+m[1]])
        x.append(iv[9+m[0]+m[1]])
        m.append(iv[10+m[0]+m[1]])
        j.append(iv[11+m[0]+m[1]:11+m[0]+m[1]+m[2]])

        lfsr1 = Generator._lfsr(c, p[0], w, m[0], j[0], x[0])
        lfsr2 = Generator._lfsr(c, p[1], w, m[1], j[1], x[1])
        lfsr3 = Generator._lfsr(c, p[2], w, m[2], j[2], x[2])

        for i in range(c):
            v = (lfsr1[i] & lfsr2[i]) ^\
                (lfsr2[i] & lfsr3[i]) ^ lfsr3[i]
            values.append(v)

        return values

    @staticmethod
    def mersenne_twister(nc, iv):
        values = []
        n = nc
        p = iv[0]
        xs = iv[-1:-p-1:-1][::-1]
        w, r, q, a, u, s, t, l, b, c = iv[1:11]

        t1 = [1]*(w-r) + [0]*r
        t2 = [0]*(w-r) + [1]*r

        d1 = int(''.join(str(x) for x in t1), 2)
        d2 = int(''.join(str(x) for x in t2), 2)

        for i in range(n):
            t12 = xs[i] & d1
            t13 = (xs[i + 1]) & d2
            Y = t12 | t13

            x = 0
            if (Y % 2 != 0):
                x = (xs[i + q] % 2 ** w) ^ (Y >> 1) ^ a
            else:
                x = (xs[i + q] % 2 ** w) ^ (Y >> 1) ^ 0

            Y = x
            Y = (Y ^ (Y >> u))
            Y = Y ^ ((Y << s) & b)
            Y = Y ^ ((Y << t) & c)
            Z = (Y ^ (Y >> l))

            xs.append(x)
            values.append(Z)
        return values


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-g', choices=['lc', 'add', '5p', 'lfsr', 'nfsr', 'mt', 'rc4', 'rsa', 'bbs'], required=True)
    parser.add_argument('-f', default='rnd.dat')
    parser.add_argument('-n', type=int, default=10000)
    parser.add_argument('-i', type=int, nargs='+', required=True)
    args = parser.parse_args()
    match args.g:
        case 'lc':
            values = Generator.linear_congruential(args.n, args.i)
            with open(args.f, 'w') as f:
                f.write('\n'.join(str(x) for x in values))
        case 'add':
            values = Generator.additive(args.n, args.i)
            with open(args.f, 'w') as f:
                f.write('\n'.join(str(x) for x in values))
        case 'lfsr':
            values = Generator.lfsr(args.n, args.i)
            with open(args.f, 'w') as f:
                f.write('\n'.join(str(x) for x in values))
        case '5p':
            values = Generator.five_parameters(args.n, args.i)
            with open(args.f, 'w') as f:
                f.write('\n'.join(str(x) for x in values))
        case 'rsa':
            values = Generator.rsa(args.n, args.i)
            with open(args.f, 'w') as f:
                f.write('\n'.join(str(x) for x in values))
        case 'bbs':
            values = Generator.bbs(args.n, args.i)
            with open(args.f, 'w') as f:
                f.write('\n'.join(str(x) for x in values))
        case 'rc4':
            values = Generator.rc4(args.n, args.i)
            with open(args.f, 'w') as f:
                f.write('\n'.join(str(x) for x in values))
        case 'nfsr':
            values = Generator.nfsr(args.n, args.i)
            with open(args.f, 'w') as f:
                f.write('\n'.join(str(x) for x in values))
        case 'mt':
            values = Generator.mersenne_twister(args.n, args.i)
            with open(args.f, 'w') as f:
                f.write('\n'.join(str(x) for x in values))


if __name__ == '__main__':
    main()
