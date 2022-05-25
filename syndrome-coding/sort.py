#!/usr/bin/env python3

import csv


codes = [
    "bch_syndrome_256_128"
]

for i in range(len(codes)):
    for t in ['ber','bler']:
        result = []
        with open(f'results/{t}_{codes[i]}.csv', 'r', newline='') as csvfile:
            reader = csv.reader(csvfile, delimiter=',')
            header = ['snr',t]
            for row in reader:
                if row[0] == 'snr':
                    continue
                result.append([float(row[0]), row[1]])
            result.sort()
            result = [header] + result
        with open(f'results/{t}_{codes[i]}.csv', 'w') as out:
            writer = csv.writer(out, delimiter=',')
            writer.writerows(result)

