import scipy
import sys
import math
import random
import numpy as np
import scipy.stats
import scipy.special
from sympy.functions.combinatorial.numbers import stirling
import matplotlib.pyplot as plt


sqrt_two = math.sqrt(2)
crit_unique_values = 49
eps = 0.001


def compute_sd(seq, l_seq, expected=False):
    if type(expected) == type(True):
        expected = compute_expected(seq, l_seq)
    return math.sqrt(sum([(seq[j] - expected) ** 2 for j in range(l_seq)]) / l_seq)


def compute_expected(seq, l_seq):
    return sum(seq) / l_seq


def uniform(seq, l_seq, d=64, targeted=True):

    seq = list(map(lambda elem : math.floor(elem * d), seq))
    alpha = .05
    critical_value = round(scipy.stats.chi2.ppf(1 - alpha, d - 1), 3)
    est, chi = l_seq / d, 0

    if targeted:
        print(f"Уровень значимости: {alpha}.")
        print(f"Критическое значение для критерия равномерности с {d - 1} степенями свободы: {critical_value}.")
        for i in range(d):
            chi += ((seq.count(i) - est) ** 2) / est
        chi = round(chi, 3)
        print(f"Вычисленное значение хи-квадрат: {chi}.")
        print(f"Представленная последовательность удовлетворяет критерию равномерности: {chi < critical_value}.")

    return seq


def series(seq, l_seq):

    d, matches = 4, {}
    d_sq = d * d
    l_unique_values = len(np.unique(seq))

    for i in range(d):
        for j in range(d):
            matches[(i, j)] = 0
    seq = uniform(seq, l_seq, d, targeted=False)
    alpha, chi = .05, 0
    critical_value, est = round(scipy.stats.chi2.ppf(1 - alpha, d_sq), 3), round(l_seq / (2 * d_sq), 3)

    for i in range(l_seq // 2):
        matches[(seq[2 * i], seq[2 * i + 1])] += 1

    print(f"Совпадения пар последовательных чисел: {matches}.")
    print(f"Теоретическое количество попаданий в каждую категорию: {est}.")
    print(f"Уровень значимости: {alpha}.")
    print(f"Параметр d принимает значение: {d}.")
    print(f"Критическое значение для критерия серий с {d_sq} степенями свободы: {critical_value}.")
    for key in matches.keys():
        chi += ((matches[key] - est) ** 2) / est
    chi = round(chi, 3)
    print(f"Вычисленное значение хи-квадрат: {chi}.")
    print(f"Представленная последовательность удовлетворяет критерию серий: {chi < critical_value}.")
    if l_unique_values < crit_unique_values:
        print(f"Из-за маленького периода ({l_unique_values}) предложенная последовательность не может удовлетворять критерию серий.")


def intervals(seq, l_seq):

    t_lower, t_upper, n_lower, n_upper = 8, 12, 1000, 2000
    _alpha = .05

    alpha, beta = random.uniform(0, 1), random.uniform(0, 1)
    while not (0.2 < abs(alpha - beta) < 0.4):
        alpha, beta = random.uniform(0, 1), random.uniform(0, 1)
    if beta < alpha: alpha, beta = beta, alpha
    alpha, beta = round(alpha, 3), round(beta, 3)
    p = beta - alpha

    s, count, chi, best_crit = 0, {}, math.inf, math.inf
    best_len, best_ints = math.inf, math. inf
    best_probs, best_count = list(), list()

    for t in range(t_lower, t_upper):
        for n in range(n_lower, n_upper):
            s, count = 0, {}
            for k in range(t + 1): count[k] = 0
            r = 0
            for j in range(l_seq):
                if alpha <= seq[j] < beta:
                    if r >= t: count[t] += 1
                    else: count[r] += 1
                    s += 1
                    if s < n: r = 0
                    else: break
                else: r += 1

            probs, chi_acc = [round(p * ((1 - p) ** _r), 3) for _r in range(t)], 0
            probs.append(round((1 - p) ** t, 4))
            critical_value = round(scipy.stats.chi2.ppf(1 - _alpha, t + 1), 3)
            for j in range(len(probs)):
                est = n * probs[j]
                if est == 0: break
                chi_acc += ((count[j] - est) ** 2) / est
            if chi_acc < chi:
                chi, best_crit = chi_acc, critical_value
                best_len, best_ints = t, n
                best_probs, best_count = probs, count

    print(f"Уровень значимости: {_alpha}.")
    print(f"Критическое значение для критерия хи-квадрат с {best_len + 1} степенями свободы: {best_crit}.")
    print(f"Найденное оптимальное значение хи-квадрата: {round(chi, 3)}.")
    print(f"Это значение было найдено при параметрах t (макс. длина интервала) = {best_len} и n (количество интервалов) = {best_ints}.")
    print(f"Используемые значения границ alpha = {alpha}, beta = {beta}.")
    print(f"Представленная последовательность удовлетворяет критерию интервалов: {chi < best_crit}.")
    print(f"Пересчитанные вероятности p_r и p_t: {best_probs}.")
    print(f"Подсчитанные значения интервалов длиной 0, 1, ..., t - 1 и >= t: {best_count}.")


def partitions(seq, l_seq, d=8, k=5):

    seq = uniform(seq, l_seq, d, targeted=False)
    groups_quan, alpha, chi, num_groups = {}, .05, 0, l_seq // k
    critical_value = round(scipy.stats.chi2.ppf(1 - alpha, k - 1), 3)
    for i in range(k, 0, -1): groups_quan[i] = 0
    counter = 0

    for i in range(num_groups):
        current_group = seq[(i * k):((i + 1) * k)]
        groups_quan[len(list(np.unique(current_group)))] += 1
    probs = [round((math.prod([i for i in range(d, d - r, -1)]) / (d ** k)) * stirling(k, r, kind=2), 4) for r in groups_quan.keys()]
    chi = sum([((value - probs[i] * num_groups) ** 2) / (probs[i] * num_groups) for i, value in enumerate(groups_quan.values())])

    print(f"Уровень значимости: {alpha}.")
    print(f"Получившийся словарь комбинаций (первый элемент соответствует количеству уникальных символов в группе): {groups_quan}.")
    print(f"Вычисленные теоретические вероятности частоты комбинаций: {probs}.")
    print(f"Критическое значение для критерия хи-квадрат с {k - 1} степенями свободы: {critical_value}.")
    print(f"Найденное значение хи-квадрат: {chi}.")
    print(f"Представленная последовательность удовлетворяет критерию разбиений: {chi < critical_value}.")


def permutations(seq, l_seq, t=4):

    df = math.prod([j for j in range(1, t + 1)])
    categ_quan = {i:0 for i in range(df)}
    alpha, num_groups = .05, l_seq // t
    critical_value = round(scipy.stats.chi2.ppf(1 - alpha, df), 3)

    def perm_analysis(perm):
        r, f = t, 0
        while r > 1:
            s = perm.index(max(perm)) + 1
            f = r * f + s - 1
            perm[r - 1], perm[s - 1] = perm[s - 1], perm[r - 1]
            perm, r = perm[:-1], r - 1

        return f

    for i in range(num_groups):
        current_group = seq[(i * t):((i + 1) * t)]
        if len(np.unique(current_group)) != t: continue
        else: categ_quan[perm_analysis(current_group)] += 1

    est, chi = round(1 / df, 3), 0
    est_cat = round(est * num_groups, 3)
    print(f"Уровень значимости: {alpha}.")
    print(f"Эмпирическое распределение по частотам: {categ_quan}.")
    print(f"Теоретическое значение вероятности для каждой категории: {est}. Теоретическое количество попаданий в каждую категорию: {est_cat}.")
    print(f"Критическое значение для критерия хи-квадрат с {df} степенями свободы: {critical_value}.")
    chi = round(sum([((categ_quan[i] - est_cat) ** 2) / est_cat for i in range(df)]), 3)
    print(f"Найденное значение хи-квадрат: {chi}.")
    print(f"Представленная последовательность удовлетворяет критерию перестановок: {chi < critical_value}.")


def monotonous(seq, l_seq):

    categ_quan, est = {}, list()
    alpha, chi, i = .05, 0, 0

    current_spree = 1
    while i < l_seq - 1:
        if seq[i] < seq[i + 1]: current_spree += 1
        else:
            if current_spree not in categ_quan:
                categ_quan[current_spree] = 1
            else: categ_quan[current_spree] += 1
            current_spree, i = 1, i + 1
        i += 1
    sum_series = sum([value for value in categ_quan.values()])
    num_categs = len(categ_quan)
    critical_value = round(scipy.stats.chi2.ppf(1 - alpha, num_categs), 3)

    fact_acc = 1
    for i in range(1, num_categs + 1):
        fact_acc *= i
        est.append(round(1 / fact_acc - 1 / (fact_acc * (i + 1)), 4))
    est = [round(sum_series * est[j], 3) for j in range(num_categs)]
    chi = round(sum([((categ_quan[j + 1] - est[j]) ** 2) / est[j] for j in range(num_categs)]), 3)

    print(f"Уровень значимости: {alpha}.")
    print(f"Эмпирические значения длин серий: {categ_quan}.")
    print(f"Теоретические значения длин серий: {est}.")
    print(f"Критическое значение для критерия хи-квадрат с {num_categs} степенями свободы: {critical_value}.")
    print(f"Найденное значение хи-квадрат: {chi}.")
    print(f"Представленная последовательность удовлетворяет критерию монотонности: {chi < critical_value}.")


def conflicts(seq, l_seq):

    s_lower, s_upper, d_lower, d_upper, n_to_print = 8, 20, 2, 8, 5
    m_lower, m_upper = 16, 128
    eps = 1e-20
    T_table = (.01, .05, .25, .50, .75, .95, .99, 1.)

    def percent_points(m, n):
        A, conflicts_et_probs = [0] * (n + 1), list()
        A[1] = j_0 = j_1 = 1

        for i in range(n - 1):
            j_1 += 1
            for j in range(j_1, j_0 - 1, -1):
                j_by_m = j / m
                A[j] = j_by_m * A[j] + (1 + 1 / m - j_by_m) * A[j - 1]
                if A[j] < eps:
                    A[j] = 0
                    if j == j_1: j_1 -= 1; continue
                    if j == j_0: j_0 += 1

        p, t, j = 0, 0, j_0 - 1
        while t != len(T_table) - 1:
            while p <= T_table[t]:
                j += 1
                p += A[j]
            conflicts_et_probs.append((n - j - 1, round(1 - p, 3)))
            t += 1

        return conflicts_et_probs

    suitables = list()
    for vec_size in range(s_lower, s_upper + 1):
        n_param = l_seq // vec_size
        for d_param in range(d_lower, d_upper + 1):
            seq_normed = uniform(seq, l_seq, d=d_param, targeted=False)
            words, n_conflicts = list(), 0
            for j_index in range(n_param):
                _slice = seq_normed[(j_index * vec_size):((j_index + 1) * vec_size)]
                if _slice not in words: words.append(_slice)
                else: n_conflicts += 1
            for m_param in range(m_lower, m_upper + 1):
                m_value = n_param * m_param
                confs_et_probs = percent_points(m_value, n_param)[::-1]
                if n_conflicts == 0 or confs_et_probs[0][0] == -1 or confs_et_probs[0][0] == 0:
                    continue
                if confs_et_probs[2][0] <= n_conflicts <= confs_et_probs[-3][0]:
                    suitables.append((confs_et_probs, n_conflicts, vec_size, n_param, d_param, m_param, m_value))

    l_suits, rand_idxs = len(suitables), list()
    for j in range(n_to_print):
        rand_idx = random.randint(0, l_suits - 1)
        if rand_idx not in rand_idxs: rand_idxs.append(rand_idx)
        else:
            while rand_idx in rand_idxs: rand_idx = random.randint(0, l_suits - 1)

    print(f"Было найдено {l_suits} наборов параметров, при которых представленная последовательность удовлетворяет критерию конфликтов.")
    print(f"Случайные {len(rand_idxs)} из них будут выведены на экран.")
    for j in range(len(rand_idxs)):
        confs_et_probs, n_conflicts, vec_size, n_param, d_param, m_param, m_value = suitables[rand_idxs[j] % l_suits]
        print(f"Размерность вектора V_j: {vec_size}. Количество векторов: {n_param}. Множитель для m: {m_param}. Само значение m: {m_value}.")
        print(f"Параметр нормирования d: {d_param}. Количество возникших конфликтов: {n_conflicts}.")
        print(f"Таблица процентных точек: {confs_et_probs}.", end='\n\n')


def chi_wrapper(seq, l_seq, intervals):

    print(f"Используемые диапазоны интервалов: {[(round(interval[0], 4), round(interval[1], 4)) for interval in intervals]}")
    ints_len = len(intervals)

    def printer(est, chi_val, chi_crit):
        print(f"Ожидаемое распределение чисел по интервалам: {est}.")
        print(f"Значение хи-квадрат: {chi_val}.")
        print(f"Результат о принятии гипотезы: {'принята' if 0 < chi_val < chi_crit else 'не принята'}.", end='\n')

    actual = [0 for i in range(ints_len)]
    for i in range(l_seq):
        for j in range(ints_len):
            if intervals[j][0] <= seq[i] < intervals[j][1]:
                actual[j] += 1; break

    est = round(l_seq / ints_len, 3)
    chi_st = round(sum([((actual[j] - est) ** 2) / est for j in range(ints_len)]), 4)
    chi_crit = round(scipy.stats.chi2.ppf(1-.05, df=(ints_len - 1)), 4)

    print(f"Количество степеней свободы: {ints_len - 1}. Уровень значимости: {0.05}.", end= ' ')
    print(f"Критическое значение хи-квадрат: {chi_crit}.")
    print(f"Наблюдаемое распределение чисел по интервалам: {actual}.", end='\n')
    printer(est, chi_st, chi_crit)


def draw_params(seq, l_seq):

    steps = (50, 100, 200, 500)
    exp_vals_stepped, _sd_vals_stepped = list(), list()

    for step in steps:
        exp_acc, exp_vals = 0, list()
        sd_vals = list()
        for j in range(l_seq // step):
            idx_l, idx_r = j * step, (j + 1) * step
            _slice_exp = compute_expected(seq[idx_l:idx_r], idx_r)
            exp_acc = exp_acc * idx_l / idx_r + _slice_exp
            exp_vals.append(exp_acc)
            sd_acc = compute_sd(seq[:idx_r], idx_r, exp_acc)
            sd_vals.append(sd_acc)
        exp_vals_stepped.append(exp_vals)
        _sd_vals_stepped.append(sd_vals)

    exp_fig, exp_axs = plt.subplots(nrows=2, ncols=2)
    exp_fig.suptitle('Сходимость мат. ожидания для различных значений шага')
    _sd_fig, _sd_axs = plt.subplots(nrows=2, ncols=2)
    _sd_fig.suptitle('Сходимость среднеквадратичного отклонения для различных значений шага')

    for j in range(len(steps)):
        j_bin = bin(j)[2:]
        l_j_bin = len(j_bin)
        stepped = [steps[j] * (k + 1) for k in range(l_seq // steps[j])]
        idx_0 = 0 if l_j_bin == 1 else int(j_bin[0])
        idx_1 = int(j_bin[0]) if l_j_bin == 1 else int(j_bin[1])

        exp_axs[idx_0, idx_1].plot(stepped, exp_vals_stepped[j])
        exp_axs[idx_0, idx_1].set_title(f"Шаг: {steps[j]}")
        _sd_axs[idx_0, idx_1].plot(stepped, _sd_vals_stepped[j])
        _sd_axs[idx_0, idx_1].set_title(f"Шаг: {steps[j]}")

    plt.show()


def main():

    file = sys.argv[1][3:]
    criterion = sys.argv[2][3:]

    f = open(file, 'r')
    seq = list(map(float, (f.read()).split(',')))
    l_seq = len(seq)
    unique_values = np.unique(seq)

    if len(unique_values) > crit_unique_values:
        seq = [seq[i] if seq[i] != 1.0 else round(np.random.uniform(), 3) for i in range(l_seq)]
        l_bord, r_bord = min(seq), max(seq)
        num_intervals = math.ceil(1 + 1.4 * math.log(l_seq))
        step = round((r_bord - l_bord) / num_intervals, 4)
        _intervals = [(l_bord + i * step, l_bord + (i + 1) * step) for i in range(num_intervals)]
    else:
        seq = [seq[i] if seq[i] != 1.0 else seq[i] - eps for i in range(l_seq)]
        unique_values = np.unique(seq)
        _intervals = [(value - eps, value + eps) for value in unique_values]

    if criterion=='a':
        chi_wrapper(seq, l_seq, _intervals)
    if criterion=='b':
        series(seq, l_seq)
    if criterion=='c':
        intervals(seq, l_seq)
    if criterion=='d':
        partitions(seq, l_seq)
    if criterion=='e':
        permutations(seq, l_seq)
    if criterion=='f':
        monotonous(seq, l_seq)
    if criterion=='g':
        conflicts(seq, l_seq)

    exp, sd = round(compute_expected(seq, l_seq), 4), round(compute_sd(seq, l_seq), 4)
    exp_th, sd_th = .5, round(math.sqrt(1 / 12), 4)
    print()
    print(f"Мат. ожидание последовательности: {exp}. Её среднеквадратичное отклонение: {sd}.")
    print(f"Эталонные значения мат. ожидания: {exp_th}; среднеквадратичного отклонения: {sd_th}.")
    rel_err_exp, rel_err_sd = round(abs(exp - exp_th) / exp, 4), round(abs(sd - sd_th) / sd, 4)
    print(f"Относительные погрешности мат. ожидания: {rel_err_exp}; среднеквадратичного отклонения: {rel_err_sd}.")

    waiter = input()
    draw_params(seq, l_seq)


if __name__ == '__main__':
    main()

