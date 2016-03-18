import fileinput

print_module = 0
for line in fileinput.input():
    if line.strip() == '>>END_MODULE':
        print_module = 0
    if print_module == 1:
        print line.strip()
    elif '>>Insert Length Distribution' in line:
        print_module = 1  