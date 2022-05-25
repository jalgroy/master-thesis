#!/usr/bin/env python3

import csv

codes = [
    "2_1_2", 
    "2_1_4", 
    "2_1_6", 
]

for i in range(len(codes)):
    for l in [128]:
        for t in ['fer','ber']:
            result = []
            with open(f'results/{t}_{codes[i]}_l{l}.csv', 'r', newline='') as csvfile:
                reader = csv.reader(csvfile, delimiter=',')
                header = ['snr',t]
                for row in reader:
                    if row[0] == 'snr':
                        continue
                    result.append([float(row[0]), row[1]])
                result.sort()
                result = [header] + result
            with open(f'results/{t}_{codes[i]}_l{l}.csv', 'w') as out:
                writer = csv.writer(out, delimiter=',')
                writer.writerows(result)

