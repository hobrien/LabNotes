import sys
import pandas as pd
import numpy as np

Combined = pd.read_table(sys.argv[1], header=None)

for index in range(4, len(Combined.columns)):
    Combined[index] = np.where(Combined[index]=='-', 0, 1)

for class_symbol in np.unique(Combined[3]):
    subset = Combined[Combined[3] == class_symbol].ix[:,4:]
    for n in range(1, len(subset.columns), 5):
        sample = []
        for x in range(1000):
            sample.append(np.count_nonzero(subset.sample(n=n, axis=1).sum(axis=1)))
        print '\t'.join((class_symbol, str(n), str(sum(sample)/100.0)))