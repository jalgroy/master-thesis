#!/usr/bin/env python3

import csv
import math


codes = [
    "uncoded_l128",
]

bob_snr = [
    6.95,
] # Tuned for ~0.1 block error rate

bob_snr_1000 = [
    9.7
]

types = ["ber","bler"]

d_bob = 10
gamma = 3

for t in types:
    for i in range(len(codes)):
        with open(f'results/{t}_{codes[i]}.csv', 'r', newline='') as csvfile:
            reader = csv.reader(csvfile, delimiter=',')
            result = [['d',t]]
            for row in reader:
                if row[0] == 'snr':
                    continue
                snr1 = bob_snr_1000[i]
                snr2 = float(row[0])
                d_eve = 10**((snr1+10*gamma*math.log10(d_bob)-snr2)/(10*gamma))
                result.append([d_eve,row[1]])
            with open(f'results/path_loss_0.001_{t}_d{d_bob}_gamma{gamma}_{codes[i]}.csv', 'w') as out:
                writer = csv.writer(out, delimiter=',')
                writer.writerows(result)

