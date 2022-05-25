#!/usr/bin/env python3

import sys
import threading

def parse_code(code):
    for l in [128]:
        print(f'{code}, length {l}')
        block_errors = {}
        bit_errors = {}
        counts = {}
        filename = f'bch_syndrome_{code}.csv'
        with open('results/' + filename, 'r', 50000000) as f:
            f.readline()
            for line in f:
                line = line.split(',')
                snr = float(line[1])
                err = int(line[2])
                if snr not in counts:
                    counts[snr] = 0
                    bit_errors[snr] = 0
                    block_errors[snr] = 0
                # Update total counts
                counts[snr] += 1
                # Update errors
                bit_errors[snr] += err
                block_errors[snr] += 0 if err == 0 else 1
        bler_f = open(f'results/bler_bch_syndrome_{code}.csv', 'w')
        ber_f = open(f'results/ber_bch_syndrome_{code}.csv', 'w')
        bler_f.write('snr,bler\n')
        ber_f.write('snr,ber\n')
        for snr in block_errors:
            block_err_rate = block_errors[snr] / counts[snr]
            bit_err_rate = bit_errors[snr] / (counts[snr]*l)
            bler_f.write(f'{snr}, {block_err_rate}\n')
            ber_f.write(f'{snr}, {bit_err_rate}\n')
        bler_f.close()
        ber_f.close()


codes = [
    "256_128" 
]

def main():
    threads = []
    
    for code in codes:
        t = threading.Thread(target=parse_code, args=(code,))
        t.start()
        threads.append(t)
    
    for t in threads:
        t.join()

if __name__ == '__main__':
    main()

