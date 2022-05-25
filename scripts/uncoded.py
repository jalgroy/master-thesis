#!/usr/bin/env python3

import math

def Q(x):
    return 1/2 * math.erfc(x / (2**0.5))

Eb = Es = 1

for l in [16,32,64,128,256]:
    with open(f'results/ber_uncoded_l{l}.csv', 'w') as ber_f, open(f'results/bler_uncoded_l{l}.csv', 'w') as bler_f:
        ber_f.write('snr,ber\n')
        bler_f.write('snr,bler\n')
        for snr in [x / 2 for x in range(-30,30)]:
            snr_lin = 10**(snr/10.0)
            N0 = Eb / snr_lin
            P_ber = Q((2*Es/N0)**0.5)
            P_bler = 1 - (1-P_ber)**l
            ber_f.write(f'{snr},{P_ber}\n')
            bler_f.write(f'{snr},{P_bler}\n')



