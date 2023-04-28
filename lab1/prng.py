import random
import sys
import strings
import math


def get_input():

    opts = {}
    args = sys.argv[1:]

    for i in range(len(args)):
        opts[args[i][1]] = args[i][3:]

    if len(opts) == 0:
        opts = {'g': 'lc', 'i': [6075, 106, 1283, random.randint(1, 1000)], 'n': 10000, 'f': 'rnd.dat'}
        return opts

    keys = opts.keys()

    if 'h' in keys:
        opts['h'] = True

    if 'g' not in keys:
        opts['g'] = 'lc'

    if 'i' in keys:
        opts['i'] = list(map(int, opts['i'].split(',')))
    else:
        opts['i'] = [6075, 106, 1283, random.randint(1, 1000)]

    if 'n' in keys:
        opts['n'] = int(opts['n'])
    else:
        opts['n'] = 10000

    if 'f' not in keys:
        opts['f'] = 'rnd.dat'

    return opts


def save_as_file(seq_num, seq_len, path):

    with open(path, 'w') as f:
        for i in range(seq_len):
            f.write(str(seq_num[i]) + ('\n' if i == seq_len - 1 else ','))


def progress_output(perc_cur, perc_25, perc_50, perc_75):

    if perc_cur == perc_25:
        print('[STATUS] Генерация чисел: 25%')
    elif perc_cur == perc_50:
        print('[STATUS] Генерация чисел: 50%')
    elif perc_cur == perc_75:
        print('[STATUS] Генерация чисел: 75%')


def rsa_machinerie(seq_len, init_x, _pow, modulo, phi, w, is_bbs=False):

    seq_num = list()
    perc_25, perc_50, perc_75 = seq_len // 4, seq_len // 2, seq_len // 4 * 3

    print('[STATUS] Генерация чисел: 0%')

    for i in range(seq_len):
        bit_seq = list()
        for j in range(w):
            init_x = pow(init_x, _pow, modulo)
            bit_seq.append(init_x % 2)
        bit_seq.reverse()
        init_x = random.randint(2, modulo - 1)
        if is_bbs:
            while math.gcd(init_x, modulo) != 1:
                init_x = random.randint(2, modulo - 1)
            inix_x = (init_x * init_x) % modulo
        if not is_bbs:
            _pow = random.randint(2, phi - 1)
            while math.gcd(_pow, phi) != 1:
                _pow = random.randint(2, phi - 1)

        seq_num.append(sum([0 if bit_seq[k] == 0 else 2 ** k for k in range(w)]))
        progress_output(i, perc_25, perc_50, perc_75)

    print('[STATUS] Генерация чисел: 100%')

    return seq_num


def bbs(IV, seq_len, path):

    p, q, n = 127, 131, 16637
    phi = (p - 1) * (q - 1)
    x, w = IV

    if math.gcd(x, n) != 1:
        print(strings.invalid_init_x)
        while math.gcd(x, n) != 1:
            x = random.randint(2, n - 1)
    x = (x * x) % n

    seq_num = rsa_machinerie(seq_len, x, 2, n, phi, w, is_bbs=True)
    save_as_file(seq_num, seq_len, path)

#python3 prng.py /g:rsa /i:514081,99991,0,127,10
def rsa(IV, seq_len, path):

    seq_num = []
    p, q, e, init_x, w = IV
    n, phi = p * q, (p - 1) * (q - 1)

    if not (1 < e < phi) or math.gcd(e, phi) != 1:
        print(strings.invalid_e)
        e = 0
        while math.gcd(e, phi) != 1:
            e = random.randint(2, phi - 1)

    if not (1 < init_x < n):
        print(strings.invalid_init_x)
        init_x = random.randint(1, n - 1)

    seq_num = rsa_machinerie(seq_len, init_x, e, n, phi, w)
    save_as_file(seq_num, seq_len, path)


"""
stream=$(python3 -c "import random; input=''.join([str(random.randint(1, 1024)) + ('' if i == 255 else ',') for i in range(256)]); print(input)")
python3 prng.py /g:rc4 /i:$stream
"""

def rc4(IV, seq_len, path):

    seq_num = []
    s_block, key, j = [i for i in range(256)], IV, 0
    perc_25, perc_50, perc_75 = seq_len // 4, seq_len // 2, seq_len // 4 * 3

    for i in range(256):
        j = ( j + s_block[i] + key[i] ) % 256
        s_block[i], s_block[j] = s_block[j], s_block[i]

    i, j = 0, 0

    print('[STATUS] Генерация чисел: 0%')

    for k in range(seq_len):

        i = ( i + 1 ) % 256
        j = ( j + s_block[i] ) % 256
        s_block[i], s_block[j] = s_block[j], s_block[i]
        t = ( s_block[i] + s_block[j] ) % 256
        seq_num.append(s_block[t])

        progress_output(k, perc_25, perc_50, perc_75)

    print('[STATUS] Генерация чисел: 100%')

    save_as_file(seq_num, seq_len, path)


mt = [0 for i in range(624)]
def twist(mt):

    lower_mask, upper_mask, reg_len = 0x7fffffff, 0x80000000, 624

    for i in range(reg_len):
        x = (mt[i] & upper_mask) + (mt[(i + 1) % reg_len] & lower_mask)
        xA = x >> 1
        if (x % 2) != 0:
            xA ^= 0x9908b0df
        mt[i] = mt[(i + 397) % reg_len] ^ xA


def extract_number(index, mt):

    if index >= 624:
        twist(mt)
        index = 0

    y = mt[index]
    y ^= ((y >> 11) & 0xffffffff)
    y ^= ((y <<  7) & 0x9d2c5680)
    y ^= ((y << 15) & 0xefc60000)
    y ^=  (y >> 18)

    index += 1

    return index, y & 0xffffffff


def MT(IV, seq_len, path):

    md, seed = IV
    register_len, seq_num = 624, []
    index, gen_num = register_len, 1812433253
    mt[0] = seed
    perc_25, perc_50, perc_75 = seq_len // 4, seq_len // 2, seq_len // 4 * 3

    for i in range(1, register_len):
        temp = gen_num * (mt[i - 1] ^ (mt[i - 1] >> 30)) + i
        mt[i] = temp & 0xffffffff

    print('[STATUS] Генерация чисел: 0%')

    for i in range(seq_len):

        index, y = extract_number(index, mt)
        seq_num.append(y % md)
        progress_output(i, perc_25, perc_50, perc_75)

    print('[STATUS] Генерация чисел: 100%')

    save_as_file(seq_num, seq_len, path)


def bit_array_conversion(input_register):

    input_register = list(map(int, str(input_register)))
    input_register.reverse()

    return input_register


def bitwise_or(fb_rev, sb_rev):

    fb_rev_s, sb_rev_s = len(fb_rev), len(sb_rev)
    min_size = min(fb_rev_s, sb_rev_s)

    return [1 if fb_rev[i] == 1 or sb_rev[i] == 1 else 0 for i in range(min_size)] + (fb_rev[min_size:] if fb_rev_s > sb_rev_s else sb_rev[min_size:])


def xor(fb_rev, sb_rev):

    fb_rev_s, sb_rev_s = len(fb_rev), len(sb_rev)
    min_size = min(fb_rev_s, sb_rev_s)

    return [(fb_rev[i] + sb_rev[i]) % 2 for i in range(min_size)] + (fb_rev[min_size:] if fb_rev_s > sb_rev_s else sb_rev[min_size:])


def nfsr(IV, seq_len, path):

    R1, R2, R3 = IV[:3]
    reg_R1, reg_R2, reg_R3, w = IV[3:]
    R1, R2, R3 = list(map(bit_array_conversion, [R1, R2, R3]))
    reg_R1, reg_R2, reg_R3 = list(map(bit_array_conversion, [reg_R1, reg_R2, reg_R3]))
    R = bitwise_or( bitwise_or( xor(R1, R2), xor(R2, R3) ), R3 )
    reg_R = bitwise_or( bitwise_or( xor(reg_R1, reg_R2), xor(reg_R2, reg_R3) ), reg_R3 )

    seq_num = lfsr_machinerie(seq_len, 0, w, [i for i in range(len(R)) if R[i] == 1], reg_R)
    save_as_file(seq_num, seq_len, path)


def lfsr_machinerie(seq_len, p, w, nonzero_coeffs, init_register=[]):

    seq_num, bit_seq = [], []
    perc_25, perc_50, perc_75 = seq_len // 4 * w, seq_len // 2 * w, seq_len // 4 * 3 * w

    print('[STATUS] Генерация чисел: 0%')

    if len(init_register) == 0:
        bit_seq = [random.randint(0, 1) for k in range(p)]
        bit_seq.reverse()
    else:
        bit_seq = init_register

    for i in range(seq_len * w):

        bit_seq.append(sum(bit_seq[i + nonzero_coeffs[k]] for k in range(len(nonzero_coeffs))) % 2)
        if i % (w - 1) == 0 and i != 0:
            current_binary = bit_seq[-w:]
            binary_converted = sum([0 if current_binary[i] == 0 else 2 ** i for i in range(w)])
            seq_num.append(binary_converted)

        progress_output(i, perc_25, perc_50, perc_75)

    print('[STATUS] Генерация чисел: 100%')

    return seq_num


def lfsr(IV, seq_len, path):

    coef_vec, init_register, w = IV
    coef_vec, init_register = list(map(int, str(coef_vec))), list(map(int, str(init_register)))
    coef_vec.reverse(), init_register.reverse()
    nonzero_coeffs = [i for i in range(len(coef_vec)) if coef_vec[i] == 1]

    seq_num = lfsr_machinerie(seq_len, 0, w, nonzero_coeffs, init_register)
    save_as_file(seq_num, seq_len, path)


def five_p(IV, seq_len, path):

    p, q_1, q_2, q_3, w = IV
    nonzero_coeffs = [0, q_1, q_2, q_3]

    seq_num = lfsr_machinerie(seq_len, p, w, nonzero_coeffs)
    save_as_file(seq_num, seq_len, path)


def add(IV, seq_len, path):

    md, k_delay, j_delay, seq_num = IV[0], IV[1], IV[2], IV[3:]
    perc_25, perc_50, perc_75 = seq_len // 4, seq_len // 2, seq_len // 4 * 3

    print('[STATUS] Генерация чисел: 0%')

    for i in range(seq_len):
        seq_num.append((seq_num[j_delay - k_delay + i] + seq_num[i]) % md)
        progress_output(i, perc_25, perc_50, perc_75) 

    print('[STATUS] Генерация чисел: 100%')

    save_as_file(seq_num[-seq_len:], seq_len, path)


def lc(IV, seq_len, path):

    seq_num = []
    md, mult, inc, init_value = IV
    seq_num.append(init_value)
    acc = init_value

    perc_25, perc_50, perc_75 = seq_len // 4, seq_len // 2, seq_len // 4 * 3

    print('[STATUS] Генерация чисел: 0%')

    for i in range(1, seq_len):
        
        acc = (mult * acc + inc) % md
        seq_num.append(acc)
        progress_output(i, perc_25, perc_50, perc_75)

    print('[STATUS] Генерация чисел: 100%')

    save_as_file(seq_num, seq_len, path)


def main():

    opts = get_input()
    method, IV, seq_len, path = opts['g'], opts['i'], opts['n'], opts['f']

    if 'h' in opts.keys():
        print(strings.help)
        return

    if method=='lc':
        lc(IV, seq_len, path)
    elif method=='add':
        add(IV, seq_len, path)
    elif method=='5p':
        five_p(IV, seq_len, path)
    elif method=='lfsr':
        lfsr(IV, seq_len, path)
    elif method=='nfsr':
        nfsr(IV, seq_len, path)
    elif method=='mt':
        MT(IV, seq_len, path)
    elif method=='rc4':
        rc4(IV, seq_len, path)
    elif method=='rsa':
        rsa(IV, seq_len, path)
    elif method=='bbs':
        bbs(IV, seq_len, path)
    else:
        print('Такого генератора нет!')
        return

if __name__ == '__main__':
    main()

#python3 prng.py /g:lc /i:1223,7,11,3
#python3 prng.py /g:add /i:100,24,55,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55 /f:'rnd_add.dat'
#python3 prng.py /g:5p /i:89,20,40,69,10
#python3 prng.py /g:lfsr /i:1000010001,1010101110100011100010101010,10 
#python3 prng.py /g:nfsr /i:10101,1000011001,10101000000001,10101001,10101000111001,1010010101010000101011100001,10 /n:10000 /f:'rnd_nfsr.dat'
#python3 prng.py /g:mt /i:2048,113 /n:1000 /f:'rnd_mt.dat'

"""
stream=$(python3 -c "import random; input=''.join([str(random.randint(1, 1024)) + ('' if i == 255 else ',') for i in range(256)]); print(input)")
python3 prng.py /g:rc4 /i:$stream /n:1000 /f:'rnd_rc4.dat'
"""
#python3 prng.py /g:rc4 /i:stream
#python3 prng.py /g:rsa /i:514081,99991,0,127,10 /n:10000 /f:'rnd_rsa.dat'
#python3 prng.py /g:bbs /i:113,10 /n:10000 /f:'bbs_rsa.dat'

#/g:lc /i:7;106;1283;6075 /n:10000 /f:rnd.dat
#/g:add /i:1000;3;15;32;23;1;54;87;543;45;24;43;867;41;21;116;6;58;43 /n:10000 /f:rnd.dat
#/g:5p /i:70;10;25;41;59;65 /n:10000 /f:rnd.dat
#/g:lfsr /i:5;5;2;1;5;12 /n:10000 /f:rnd.dat
#/g:nfsr /i:2;10;5;9;3;1;2;5;9;233;4;1;5;7;9;3;10;01;00 /n:10000 /f:rnd.dat
#/g:ms /i:5;4;2;3;13;2;1;1;3;4;6;4;7;5;2;9 /n:1000 /f:rnd.dat
#/g:rc4 /i:9;86;114;217;41;207;80;50;231;180;228;211;103;48;90;105;120;237;240;253;118;187;230;21;144;117;26;252;226;87;191;19;65;71;236;147;248;36;84;245;178;208;210;64;176;42;66;166;51;7;12;222;206;162;229;18;129;9;13;34;213;101;216;20;89;112;56;215;171;95;239;72;94;30;106;141;63;33;235;138;93;212;85;104;61;238;44;10;193;70;124;116;81;68;243;99;155;74;24;196;160;125;158;165;122;205;29;223;140;186;254;249;188;137;123;174;173;185;214;246;115;46;14;17;11;60;37;164;3;97;6;35;108;119;43;153;146;255;200;136;31;58;5;27;151;57;242;22;190;109;47;224;143;184;232;227;149;163;96;110;142;82;111;49;203;69;98;126;38;113;28;15;195;55;67;40;159;88;152;100;194;102;52;250;218;1;167;169;0;83;189;179;247;177;181;168;75;198;92;91;134;148;161;241;32;183;130;77;175;78;2;133;127;39;244;221;54;234;8;139;201;251;62;4;192;16;107;157;59;121;204;156;25;145;73;202;172;154;131;128;220;79;209;219;197;132;225;45;23;135;233;53;199;76;150;170;182 /n:10000 /f:rnd.dat
#/g:rsa /i:211;4;15;5;4 /n:10000 /f:rnd.dat
#/g:bbs /i:321;50;10 /n:10000 /f:rnd.dat