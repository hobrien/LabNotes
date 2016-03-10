import sys, warnings, subprocess

line1 = ''
line2 = ''

for line in sys.stdin:
    subprocess.call(['ln', '-s', line.strip(), "/home/heath/Raw/%s.fastq.gz" % line.strip().split('/')[-1].split('.')[0]])
