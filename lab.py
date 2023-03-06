from functools import reduce
import math
import random


class CriterionChecker:

    def __init__(self):
        self.sequence = [
            random.uniform(0.0, 1.0) for _ in range(10000)
        ]
        with open('gen.txt', 'w') as f:
            f.write('\n'.join(f'{x:.4f}' for x in self.sequence))

        self.sequence_aux = [
            random.randint(1000, 100000) for _ in range(10000)
        ]

    @staticmethod
    def mean(seq):
        return reduce(lambda x, y: x+y, seq) / len(seq)

    @staticmethod
    def variance(seq):
        me = CriterionChecker.mean(seq)
        return reduce(lambda x, y: x+y,
                      map(lambda x: (x-me)**2, seq)) / (len(seq)-1)

    @staticmethod
    def standard_deviation(seq):
        return math.sqrt(CriterionChecker.variance(seq))

    @staticmethod
    def chi_squared(seq):
        n = len(seq)
        min_v = len(list(filter(lambda x: x < 0.5, seq)))
        max_v = len(list(filter(lambda x: x > 0.5, seq)))
        s1 = (min_v - n/2) ** 2
        s2 = (max_v - n/2) ** 2
        s = (s1 / (n/2)) + (s2 / (n/2))
        return s

    @staticmethod
    def series(seq):
        seq_ = sorted(seq)
        n = len(seq_)
        res_1 = 0
        res_2 = 0
        res_3 = 0
        res_4 = 0
        for i in range(len(seq_) - 1):
            val1 = seq_[i]
            val2 = seq_[i + 1]
            if ((val1 > 0) and (val1 < 1/4)) and ((val2 > 0) and (val2 < 1/4)):
                res_1 += 1
            elif ((val1 > 1/4) and (val1 < 2/4)) and ((val2 > 1/4) and (val2 < 2/4)):
                res_2 += 1
            if ((val1 > 2/4) and (val1 < 3/4)) and ((val2 > 2/4) and (val2 < 3/4)):
                res_3 += 1
            elif ((val1 > 3/4) and (val1 < 1)) and ((val2 > 3/4) and (val2 < 1)):
                res_4 += 1
        count = n / 4
        d1 = (res_1 - count) ** 2
        d2 = (res_2 - count) ** 2
        d3 = (res_3 - count) ** 2
        d4 = (res_4 - count) ** 2
        s = (d1 / (count)) + (d2 / (count)) + \
            (d3 / (count)) + (d4 / (count))
        return s

    @staticmethod
    def interval(seq):
        t = 2
        r = 0
        count = [0, 0, 0]
        for el in seq:
            if ((el < 0.2) and (el > 0)) or ((el > 0.6) and (el < 1)):
                r += 1
            if ((el > 0.2) and (el < 0.6)):
                if (r >= t):
                    count[t] += 1
                    r = 0
                else:
                    count[r] += 1
                    r = 0
        am = sum(count)
        el1 = am * 0.4
        el2 = am * 0.24
        el3 = am * 0.36
        d1 = (count[0] - el1) ** 2
        d2 = (count[1] - el2) ** 2
        d3 = (count[2] - el3) ** 2
        s = d1/el1 + d2/el2 + d3/el3
        return s

    @staticmethod
    def transform(seq):
        res = []
        diff = 0
        d_same = 0
        t_same = 0
        for i in range(len(seq)):
            p = str(seq[i])
            if (p == '0'):
                p = "100"
            if (p[0] != p[1]) and (p[0] != p[2]) and (p[1] != p[2]):
                diff += 1
            if ((p[0] == p[1]) and (p[0] != p[2])) or ((p[0] == p[2]) and (p[0] != p[1])) or ((p[1] == p[2]) and (p[0] != p[1])):
                d_same += 1
            if (p[0] == p[1]) and (p[0] == p[2]) and (p[1] == p[2]):
                t_same += 1
        res.append(diff)
        res.append(d_same)
        res.append(t_same)
        return res

    @staticmethod
    def lef_chi(x, n):
        s = 0
        am1 = 0.72 * (n)
        am2 = 0.27 * (n)
        am3 = 0.01 * (n)
        d1 = (x[0] - am1) ** 2
        d2 = (x[1] - am2) ** 2
        d3 = (x[2] - am3) ** 2
        s = (d1 / (am1)) + (d2 / (am2)) + (d3 / (am3))
        return s

    @staticmethod
    def splitting(seq):
        seq_ = []
        for el in seq:
            d = int(el * 1000)
            if(len(str(d)) < 3):
                d *= 10
            if(len(str(d)) < 3):
                d *= 10
            seq_.append(d)
        p_p = CriterionChecker.transform(seq_)
        res = CriterionChecker.lef_chi(p_p, len(seq_))
        return res

    @staticmethod
    def permutation(seq):
        n = len(seq)
        r1 = 0
        r2 = 0
        for i in range(len(seq) - 1):
            val1 = seq[i]
            val2 = seq[i + 1]
            if (val1 < 0.5) and (val2 < 0.5):
                r1 += 1
            if (val1 > 0.5) and (val2 > 0.5):
                r2 += 1
        am = n / 4
        d1 = (r1 - am) ** 2
        d2 = (r2 - am) ** 2
        s = d1/am + d2/am
        return s

    @staticmethod
    def monotone(seq):
        t = 2
        r = 0
        count = [0, 0, 0]
        for i in range(len(seq) - 1):
            if (seq[i] < seq[i + 1]):
                r += 1
            if (seq[i] > seq[i + 1]):
                i += 1
                if (r >= t + 1):
                    count[t] += 1
                    r = 1
                else:
                    count[r - 1] += 1
                    r = 1
        am = sum(count)
        am1 = am * 0.34
        am2 = am * 0.40
        am3 = am * 0.26
        d1 = (count[0] - am1) ** 2
        d2 = (count[1] - am2) ** 2
        d3 = (count[2] - am3) ** 2
        s = (d1/am1 + d2/am2 + d3/am3) / 10
        return s

    @staticmethod
    def conflict(seq):
        n = 10000
        m = 100000
        teor_conf = (n**2) / (2*m)
        pre_set = set(seq)
        conf_am = len(seq) - len(pre_set)
        am = n - teor_conf
        d1 = (conf_am - teor_conf) ** 2
        d2 = (len(pre_set) - am) ** 2
        s = (d1/teor_conf + d2/am) / 100
        return s


def main():
    cc = CriterionChecker()
    print("Математическое ожидание ", CriterionChecker.mean(cc.sequence))
    print("Дисперсия", CriterionChecker.variance(cc.sequence))
    print("Среднеквадратичное отклонение",
          CriterionChecker.standard_deviation(cc.sequence))

    chi = CriterionChecker.chi_squared(cc.sequence)
    if ((chi > 0.0158) and (chi < 3.8415)):
        print("Критерий Хи-квадрат подтвержден  " + str(chi))
    else:
        print("Критерий Хи-квадрат не подтвержден  " + str(chi))

    series = CriterionChecker.series(cc.sequence)
    if ((series > 0.5844) and (series < 7.8147)):
        print("Критерий Cерий подтвержден  " + str(series))
    else:
        print("Критерий Cерий не подтвержден  " + str(series))

    interval = CriterionChecker.interval(cc.sequence)
    if ((interval > 0.5844) and (interval < 7.8147)):
        print("Критерий Интервалов подтвержден  " + str(interval))
    else:
        print("Критерий Интервалов не подтвержден  " + str(interval))

    split = CriterionChecker.splitting(cc.sequence)
    if ((split > 0.2107) and (split < 5.9915)):
        print("Критерий Разбиений подтвержден  " + str(split))
    else:
        print("Критерий Разбиений не подтвержден  " + str(split))

    permut = CriterionChecker.permutation(cc.sequence)
    if ((permut > 0.0158) and (permut < 3.8415)):
        print("Критерий Перестановок подтвержден  " + str(permut))
    else:
        print("Критерий Перестановок не подтвержден  " + str(permut))

    monotone = CriterionChecker.monotone(cc.sequence)
    if ((monotone > 0.2107) and (monotone < 5.9915)):
        print("Критерий Монотонности подтвержден  " + str(monotone))
    else:
        print("Критерий Монотонности не подтвержден  " + str(monotone))

    confl = CriterionChecker.conflict(cc.sequence_aux)
    if ((confl > 0.0158) and (confl < 3.8415)):
        print("Критерий Кофликтов подтвержден  " + str(confl))
    else:
        print("Критерий Конфликтов не подтвержден  " + str(confl))


if __name__ == '__main__':
    main()
