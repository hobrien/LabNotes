#!/usr/bin/env python
import sys

read1 = []
read2 = []
for i in range(1, len(sys.argv)-1, 2):
    read1 += sys.argv[i:i+1]
    read2 += sys.argv[i+1:i+2]
    
print ','.join(read1), ','.join(read2)

        
