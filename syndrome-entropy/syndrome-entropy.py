#!/usr/bin/env python3

'''Implementation of equivocation calculation for syndrome codes from Al-Hassan, Ahmed, and Tomlinson 2014.
'''

import numpy as np
from numpy.polynomial import Polynomial

H0 = [1,2,4,8,16,32,64,126,128,256,512,826,1024,2048,3879,4096,7163,7913,8192,9215,9632,10552,16384,16975,17378,17779,18843,19664,21136,21973,22578,23393,24092,24495,25144,26321,26409,26640,26663,27411,28092,28622,29302,29977,31397,31871,32395,32607]
H1 = [1,2,4,8,16,32,64,128,256,512,1024,2048,3879,4096,5541,7031,7160,7913,8192,9215,13987,14289,16384,16975,17378,17579,18413,18843,18960,19350,19955,21973,23259,23393,24092,24495,25609,25698,26321,26409,27411,28092,28133,28622,28816,31397,31675,31871,32153]

H_BCH_16_32 = [1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096, 8192, 16384, 32768, 60387, 15398, 30796, 61592, 2771, 5541, 11081, 22161, 44321, 45474, 34983, 64174, 7871, 15741, 31481, 62961]
H_BCH_32_64 = [1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096, 8192, 16384, 32768, 65536, 131072, 262144, 524288, 1048576, 2097152, 4194304, 8388608, 16777216, 33554432, 67108864, 134217728, 268435456, 536870912, 1073741824, 2147483648, 2432887841, 3003593824, 4144792801, 2132534752, 4265069504, 1835891617, 3671783233, 615973536, 1231947072, 2463894144, 3032152353, 4168556128, 1642909921, 3285819841, 380743584, 761487168, 1522974336, 3045948672, 4195882529, 1697824864, 3395649728, 97119649, 194239297, 388478593, 776957185, 1553914369, 3107828737, 3816391712, 1475470433, 2950940865, 3469049248, 210245473]

c0 = (H0, 15)
c1 = (H1, 15)
c2 = (H_BCH_16_32, 16)
c3 = (H_BCH_32_64, 32)


def xor_mult(p1,c,e,m):
    res = np.zeros(2**m, dtype=np.float32)
    for i,p in enumerate(p1):
        res[i ^ e] = p*c
    return Polynomial(res)

def equivocation(H, m, p_e):
    n = len(H)
    pz_s = Polynomial([1])
    q = 1-p_e
    
    for i in range(n):
        qs = q*pz_s
        pz_s = xor_mult(pz_s, p_e, H[i],m)
        pz_s = pz_s + qs
    
    Eq = 0
    for beta in pz_s:
        if beta != 0:
            Eq += beta*np.log2(beta)
    
    return -Eq

def main():
    #eq1 = equivocation(H0, 15, 0.05)
    #print(f'Equivocation of m=15, n=48 code with p_e=0.05: {eq1}')
    #print(f'\t Equivocation rate: {eq1/15}')

    print('p,eq')
    print('0,0')
    for p in [x / 80.0 for x in range(1, 40)]:
        eq = equivocation(c2[0], c2[1], p)
        print(f'{p},{eq}')

if __name__ == '__main__':
    main()

