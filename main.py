import random
import math
import matplotlib.pyplot as plt
# from array import *
# import numpy as np
import random

n = 10000
pre = []
math_ozid = 0


def toFixed(f: float, n=0):
    a, b = str(f).split('.')
    return '{}.{}{}'.format(a, b[:n], '0'*(n-len(b)))


# print(toFixed(1234.3434, 3))


i = 0
n = 10000
file2 = open('generator.txt', 'w')
file3 = open('new_generator.txt', 'w')
while i < n:
    # uniform() - возвращает случайное вещественное число в указанном промежутке
    x = random.uniform(0.0, 1.0)
    # uniform() - возвращает случайное интовое число в опредленном промежутке
    y = random.randint(1000, 100000)
    if (y < 10000):
        y = str('0' + str(y))
    file2.write(str(toFixed(x, 3)) + '\n')
    file3.write(str(y) + '\n')
    #print(str(toFixed(x, 3)))
    # print(str(y))
    i = i + 1


def variance(arr, math, n):
    d = 0
    for i in arr:
        d = ((i - math) ** 2) + d
    d = d / (n - 1)
    return d


def gen(fn):
    f = open(fn, "r")
    rez = []
    for line in f:
        rez.append(float(line))
    f.close()
    return rez


def hi_2(x, n):
    min_v = 0
    max_v = 0
    s = 0
    for key in x:
        if (key < 0.5):
            min_v += 1
        else:
            max_v += 1
    s1 = (min_v - n/2) ** 2
    s2 = (max_v - n/2) ** 2
    s = (s1 / (n/2)) + (s2 / (n/2))
    s = round(s, 3)
    return s


def series(x):
    x = sorted(x)
    n = len(x)
    res_1 = 0
    res_2 = 0
    res_3 = 0
    res_4 = 0
    for i in range(len(x) - 1):
        val1 = x[i]
        val2 = x[i + 1]
        if ((val1 > 0) and (val1 < 1/4)) and ((val2 > 0) and (val2 < 1/4)):
            res_1 = res_1 + 1
        else:
            if ((val1 > 1/4) and (val1 < 2/4)) and ((val2 > 1/4) and (val2 < 2/4)):
                res_2 = res_2 + 1
        if ((val1 > 2/4) and (val1 < 3/4)) and ((val2 > 2/4) and (val2 < 3/4)):
            res_3 = res_3 + 1
        else:
            if ((val1 > 3/4) and (val1 < 1)) and ((val2 > 3/4) and (val2 < 1)):
                res_4 = res_4 + 1
    kol = n / 4
    tmp1 = (res_1 - kol) ** 2
    tmp2 = (res_2 - kol) ** 2
    tmp3 = (res_3 - kol) ** 2
    tmp4 = (res_4 - kol) ** 2
    s = (tmp1 / (kol)) + (tmp2 / (kol)) + (tmp3 / (kol)) + (tmp4 / (kol))
    s = round(s, 3)
    return s


def interval(x):
    a = 0.2
    b = 0.6
    t = 2
    count = []
    r = 0
    kol = 0
    for key in range(3):
        count.append(0)
    for key in x:
        if ((key < 0.2) and (key > 0)) or ((key > 0.6) and (key < 1)):
            r = r + 1
        if ((key > 0.2) and (key < 0.6)):
            if (r >= t):
                count[t] = count[t] + 1
                r = 0
            else:
                count[r] = count[r] + 1
                r = 0
    for key in count:
        kol = kol + key
    kol1 = kol * 0.4
    kol2 = kol * 0.24
    kol3 = kol * 0.36
    tmp1 = (count[0] - kol1) ** 2
    tmp2 = (count[1] - kol2) ** 2
    tmp3 = (count[2] - kol3) ** 2
    s = (tmp1 / (kol1)) + (tmp2 / (kol2)) + (tmp3 / (kol3))
    s = round(s, 3)
    return s


def preobr(x):
    rez = []
    razn_1 = 0
    d_same = 0
    t_same = 0
    for i in range(len(x)):
        p = str(x[i])
        if (p == '0'):
            p = "100"
        if (p[0] != p[1]) and (p[0] != p[2]) and (p[1] != p[2]):
            razn_1 += 1
        if ((p[0] == p[1]) and (p[0] != p[2])) or ((p[0] == p[2]) and (p[0] != p[1])) or ((p[1] == p[2]) and (p[0] != p[1])):
            d_same += 1
        if (p[0] == p[1]) and (p[0] == p[2]) and (p[1] == p[2]):
            t_same += 1
    rez.append(razn_1)
    rez.append(d_same)
    rez.append(t_same)
    return rez


def lef_chi(x, n):
    s = 0
    kol1 = 0.72 * (n)
    kol2 = 0.27 * (n)
    kol3 = 0.01 * (n)
    tmp1 = (x[0] - kol1) ** 2
    tmp2 = (x[1] - kol2) ** 2
    tmp3 = (x[2] - kol3) ** 2
    s = (tmp1 / (kol1)) + (tmp2 / (kol2)) + (tmp3 / (kol3))
    return s


def splitting(x):
    new_x = []
    for key in x:
        tmp = int(key * 1000)
        if(len(str(tmp)) < 3):
            tmp = tmp * 10
        if(len(str(tmp)) < 3):
            tmp = tmp * 10
        new_x.append(tmp)
    p_p = preobr(new_x)
    rez = lef_chi(p_p, len(new_x))
    rez = round(rez, 3)
    return rez


def permutation(x):
    n = len(x)
    res_1 = 0
    res_2 = 0
    for i in range(len(x) - 1):
        val1 = x[i]
        val2 = x[i + 1]
        if (val1 < 0.5) and (val2 < 0.5):
            res_1 += 1
        if (val1 > 0.5) and (val2 > 0.5):
            res_2 += 1
        i = i + 1
    kol1 = n / 4
    tmp1 = (res_1 - kol1) ** 2
    tmp2 = (res_2 - kol1) ** 2
    s = (tmp1 / (kol1)) + (tmp2 / (kol1))
    s = round(s, 3)
    return s


def monotonnost(x):
    t = 2
    count = []
    r = 0
    kol = 0
    for key in range(3):
        count.append(0)
    for key in range(len(x) - 1):
        if (x[key] < x[key + 1]):
            r = r + 1
        if (x[key] > x[key + 1]):
            key = key + 1
            if (r >= t + 1):
                count[t] = count[t] + 1
                r = 1
            else:
                count[r - 1] = count[r - 1] + 1
                r = 1
    for key in count:
        kol = kol + key
    kol1 = kol * 0.34
    kol2 = kol * 0.40
    kol3 = kol * 0.26
    tmp1 = (count[0] - kol1) ** 2
    tmp2 = (count[1] - kol2) ** 2
    tmp3 = (count[2] - kol3) ** 2
    s = ((tmp1 / (kol1)) + (tmp2 / (kol2)) + (tmp3 / (kol3))) / 10
    s = round(s, 3)
    return s


def conflict():
    n = 10000
    m = 100000
    f = open("new_generator.txt", "r")
    x = []
    for line in f:
        x.append(int(line))
    f.close()
    kol_vo_conf = 0
    teor_conf = (n ** 2) / (2 * m)
    pre_set = set(x)
    kol_vo_conf = len(x) - len(pre_set)
    kol1 = n - teor_conf
    tmp1 = (kol_vo_conf - teor_conf) ** 2
    tmp2 = (len(pre_set) - kol1) ** 2
    s = ((tmp1 / teor_conf) + (tmp2 / kol1)) / (100)
    return(s)


def statistic():
    n = 1000
    math_ozid = 0
    x = []
    arr_m = []
    arr_q = []
    arr_n = []
    for i in range(100):
        for key in range(n):
            a = random.random()
            a = round(a, 3)
            x.append(a)
        for key in x:
            math_ozid += key
        math_ozid = math_ozid / len(x)
        math_ozid = round(math_ozid, 3)
        q_otklon = math.sqrt(variance(x, math_ozid, len(x)))
        q_otklon = round(q_otklon, 3)
        arr_m.append(math_ozid)
        arr_q.append(q_otklon)
        arr_n.append(n)
        n += 1000
        math_ozid = 0
        x = []
    m1 = arr_m[0]
    q1 = arr_q[0]
    print("Мат.ожидание = " + str(m1))
    print("Среднеквадратичное отклонение = " + str(q1))
    pogreshnost_m = (math.fabs(m1 - 0.5)) / m1
    pogreshnost_q = (math.fabs(q1 - (math.sqrt(1 / 12)))) / q1
    pogreshnost_m = round(pogreshnost_m, 3)
    pogreshnost_q = round(pogreshnost_q, 3)
    #print("Погрешность для мат.ожидания в % = " + str(100 * pogreshnost_m))
    #print("Погрешность для среднеквадратичного отклонения в % = " + str(100 * pogreshnost_q))
    print(' ')
    fig = plt.figure()
    for i in range(len(arr_m)):
        plt.scatter(arr_m[i], arr_n[i], c='b', s=1)
    plt.show()
    for i in range(len(arr_q)):
        plt.scatter(arr_q[i], arr_n[i], c='b', s=1)
    plt.show()


print('preobr', preobr([0.1, 0.2, 0.3, 0.4]))

statistic()

x = gen("generator.txt")

# Мат.ожидание
for key in x:
    math_ozid += key
math_ozid = math_ozid / len(x)
math_ozid = round(math_ozid, 3)
#print("Мат.ожидание = " + str(math_ozid))

# Среднеквадратичное отклонение
q_otklon = math.sqrt(variance(x, math_ozid, len(x)))
q_otklon = round(q_otklon, 3)
print("Среднеквадратичное отклонение = " + str(q_otklon))

# Критерий Хи-Квадрат
hi = hi_2(x, len(x))
if ((hi > 0.0158) and (hi < 3.8415)):
    print("Критерий Хи-квадрат подтвержден и он равен " + str(hi))
else:
    print("Критерий Хи-квадрат не подтвержден и он равен " + str(hi))

# Критерий Серий
series = series(x)
if ((series > 0.5844) and (series < 7.8147)):
    print("Критерий Cерий подтвержден и оно равен " + str(series))
else:
    print("Критерий Cерий не подтвержден и оно равен " + str(series))

# Критерий Интервалов
interval = interval(x)
if ((interval > 0.5844) and (interval < 7.8147)):
    print("Критерий Интервалов подтвержден и он равен " + str(interval))
else:
    print("Критерий Интервалов не подтвержден и он равен " + str(interval))

# Критерий Разбиений
razbien = splitting(x)
if ((razbien > 0.2107) and (razbien < 5.9915)):
    print("Критерий Разбиений подтвержден и он равен " + str(razbien))
else:
    print("Критерий Разбиений не подтвержден и он равен " + str(razbien))

# Критерий Перестановок
permut = permutation(x)
if ((permut > 0.0158) and (permut < 3.8415)):
    print("Критерий Перестановок подтвержден и он равен " + str(permut))
else:
    print("Критерий Перестановок не подтвержден и он равен " + str(permut))

# Критерий Монотонности
monoton = monotonnost(x)
if ((monoton > 0.2107) and (monoton < 5.9915)):
    print("Критерий Монотонности подтвержден и он равен " + str(monoton))
else:
    print("Критерий Монотонности не подтвержден и он равен " + str(monoton))

# Критерий Конфликтов
confl = conflict()
if ((confl > 0.0158) and (confl < 3.8415)):
    print("Критерий Кофликтов подтвержден и он равен " + str(confl))
else:
    print("Критерий Конфликтов не подтвержден и он равен " + str(confl))
