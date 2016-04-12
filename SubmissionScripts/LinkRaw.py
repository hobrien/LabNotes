import sys, warnings, subprocess

line1 = ''
line2 = ''

for line in sys.stdin:
    #subprocess.call(['ln', '-s', line.strip(), "/home/mpmho/Raw/%s.fastq.gz" % line.strip().split('/')[-1].split('.')[0]])
    subprocess.call(['ln', '-s', line.strip(), "/home/mpmho/Temp/%s.fastq" % line.strip().split('/')[-1].split('.')[0]])
