import argparse
from functools import reduce
import math


class ReadWrite():

    def __init__(self, input_file, output_file, delimiter='\n'):
        self.input_file = input_file
        self.output_file = output_file
        self.delimiter = delimiter

    def read_values_from_file(self):
        with open(self.input_file, 'r') as f:
            values = list(map(int, f.read().split(self.delimiter)))
        return values

    def write_values_to_file(self, values):
        with open(self.output_file, 'w') as f:
            f.write(f'{self.delimiter}'.join(f'{x:.4f}' for x in values))


class Distributor():

    @staticmethod
    def st(x, a, b, m):
        return [b*(x[i]/m) + a
                for i in range(len(x))]

    @staticmethod
    def tr(x, a, b, m):
        return [a + b*((x[i-1]/m) + (x[i]/m) - 1)
                for i in range(1, len(x))]

    @staticmethod
    def ex(x, a, b, m):
        return [a - b*math.log(x[i]/m)
                for i in range(len(x))]

    @staticmethod
    def nr(x, a, b, m):
        values = []
        for i in range(0, len(x), 2):
            y1 = a + b * math.sqrt(
                math.fabs(-2 * math.log(
                    math.fabs(1 - x[i]/m))
                )) * math.cos(2 * math.pi * x[i+1]/m)
            y2 = a + b * math.sqrt(
                math.fabs(-2 * math.log(
                    math.fabs(1 - x[i]/m))
                )) * math.sin(2 * math.pi * x[i+1]/m)
            values.append(y1)
            values.append(y2)

        return values

    @staticmethod
    def ln(x, a, b, m):
        return [a + math.exp(b - x[i]/m)
                for i in range(len(x))]

    @staticmethod
    def ls(x, a, b, m):
        return [a + b*math.log(math.fabs((x[i]/m) / (1 - x[i]/m)))
                for i in range(len(x))]

    @staticmethod
    def _bin_coefficient(b, k):
        return math.factorial(b) /\
            (math.factorial(k) * math.factorial(b - k))

    @staticmethod
    def bi(x, a, b, m):
        b = int(b)
        values = []
        for i in range(len(x)):
            u = x[i] / m
            y = 0
            s = 0
            k = 0

            while (True):
                s = s + \
                    Distributor._bin_coefficient(
                        b, k) * (a**k) * ((1-a) ** (b-k))
                if s > u:
                    y = k
                    break
                if k < b - 1:
                    k = k + 1
                    continue
                y = b
                break
            values.append(y)

        return values

    @staticmethod
    def gm(x, a, b, m, c):
        values = []
        if c - int(c) == 0.5:
            k = int(c-0.5)
            z = Distributor.nr(x, 0, 1, m)
            for i in range(0, len(x), k):
                sliced = x[i:i+k]
                lg = math.log(
                    reduce(lambda x, y: x*y,
                           map(lambda x: math.fabs(1 - x/m),
                               sliced)))
                y = a + b*((z[i]**2)//2 - lg)
                values.append(y)
        else:
            for i in range(0, len(x), int(c)):
                sliced = x[i:i+int(c)]
                lg = math.log(
                    reduce(lambda x, y: x*y,
                           map(lambda x: math.fabs(1 - x/m),
                               sliced)))
                y = a - b*lg
                values.append(y)

        return values


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-d', choices=['st', 'ex', 'tr', 'nr', 'gm', 'ln', 'ls', 'bi'], required=True)
    parser.add_argument('-f', default='input.dat')
    parser.add_argument('-o', default='output.dat')
    parser.add_argument('-a', type=float, required=True)
    parser.add_argument('-b', type=float, required=True)
    parser.add_argument('-c', type=float)
    parser.add_argument('-m', type=int, required=True)
    args = parser.parse_args()
    rw = ReadWrite(args.f, args.o)
    match args.d:
        case 'st':
            x = rw.read_values_from_file()
            values = Distributor.st(x, args.a, args.b, args.m)
            rw.write_values_to_file(values)
        case 'tr':
            x = rw.read_values_from_file()
            values = Distributor.tr(x, args.a, args.b, args.m)
            rw.write_values_to_file(values)
        case 'ex':
            x = rw.read_values_from_file()
            values = Distributor.ex(x, args.a, args.b, args.m)
            rw.write_values_to_file(values)
        case 'nr':
            x = rw.read_values_from_file()
            values = Distributor.nr(x, args.a, args.b, args.m)
            rw.write_values_to_file(values)
        case 'ln':
            x = rw.read_values_from_file()
            values = Distributor.ln(x, args.a, args.b, args.m)
            rw.write_values_to_file(values)
        case 'ls':
            x = rw.read_values_from_file()
            values = Distributor.ls(x, args.a, args.b, args.m)
            rw.write_values_to_file(values)
        case 'bi':
            x = rw.read_values_from_file()
            values = Distributor.bi(x, args.a, args.b, args.m)
            rw.write_values_to_file(values)
        case 'gm':
            x = rw.read_values_from_file()
            values = Distributor.gm(x, args.a, args.b, args.m, args.c)
            rw.write_values_to_file(values)


if __name__ == '__main__':
    main()
