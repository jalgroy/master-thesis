#!/usr/bin/env python3

import sys

def main():
    filename = sys.argv[1]
    mode = 0
    if len(sys.argv) > 2:
        mode = 0 if sys.argv[2] == 'ber' else 1

    with open(filename, 'r', 50000000) as f:
        f.readline()
        if mode == 0:
            print('snr,ber')
        else:
            print('snr,fer')
        for line in f:
            line = line.split(',')
            l = int(line[0])
            snr = float(line[1])
            n = int(line[2])
            block_errors = int(line[3])
            bit_errors = int(line[4])
            if mode != 0 and (block_errors < 100 or n-block_errors < 100):
                continue
            fer = block_errors / n
            ber = bit_errors / (n*l)
            if mode == 0:
                print(f'{snr},{ber}')
            else:
                print(f'{snr},{fer}')

if __name__ == '__main__':
    main()

