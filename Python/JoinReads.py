#!/usr/bin/env python
import sys

read1 = ''
read2 = ''
for i in range(1, len(sys.argv)-1, 2):
    read1 = read1 + ',' + sys.argv[i]
    read2 = read2 + ',' + sys.argv[i+1]
    
print read1, read2

        
