#!/usr/bin/env python3

import math

def h(p):
    return -p*math.log2(p)-(1-p)*math.log2(1-p)

def Q(x):
    return 1/2 * math.erfc(x / (2**0.5))

Eb = Es = 1

for l in [16,32,64,128,256]:
    with open(f'results/H_uncoded_l{l}.csv', 'w') as H_f:
        H_f.write('snr,H,d\n')
        for snr in [x / 2 for x in range(-40,40)]:
            snr_lin = 10**(snr/10.0)
            N0 = Eb / snr_lin
            P_ber = Q((2*Es/N0)**0.5)
            H = l*h(P_ber)
            H_f.write(f'{snr},{H},{H/l}\n')
