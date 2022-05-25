#!/usr/bin/env python3

import csv

codes = [
    #"32_16",
    #"64_32",
    #"128_64"
    "256_128"
]

bob_snr_10 = [
    8.68,
    9.35,
    9.95
] # Tuned for ~0.1 block error rate

bob_snr_100 = [
    11.93
]

bob_snr_1000 = [
    13
]

for i in range(len(codes)):
    with open(f'results/bler_bch_syndrome_{codes[i]}.csv', 'r', newline='') as csvfile:
        reader = csv.reader(csvfile, delimiter=',')
        result = [['snr','bler']]
        for row in reader:
            if row[0] == 'snr':
                continue
            result.append([str(bob_snr_1000[i]-float(row[0])), row[1]])
        with open(f'results/tuned_0.001_bler_bch_syndrome_{codes[i]}.csv', 'w') as out:
            writer = csv.writer(out, delimiter=',')
            writer.writerows(result)

