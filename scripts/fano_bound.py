#!/usr/bin/env python3

import sys
import csv

from math import log2

def h(p):
    return -p*log2(p)-(1-p)*log2(1-p)

def main():
    if len(sys.argv) != 3:
        print('Usage: ./fano_bound.py <filename> <bit_len>')
        sys.exit()

    filename = sys.argv[1]
    bit_len = int(sys.argv[2])
    with open(filename, 'r') as f:
        reader = csv.reader(f, delimiter=',')
        for row in reader:
            if row[0] == 'snr':
                print('snr,Hmax')
                continue
            snr = row[0]
            Pe = float(row[1])
            if Pe == 0.0:
                H = 0
            elif Pe == 1.0:
                H = Pe*log2(2**bit_len - 1)
            else:
                H = h(Pe) + Pe*log2(2**bit_len - 1)
            print(f'{snr},{H}')

if __name__ == '__main__':
    main()

