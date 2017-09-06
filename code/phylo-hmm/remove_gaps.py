import sys

x = sys.argv[1].find('.')
prefix = sys.argv[1][:x]
suffix = sys.argv[1][x:]
f = open(sys.argv[1], 'r')
fout = open(prefix + '_nogaps' + suffix, 'w')
line = f.readline()
headers = []
seqs = []
while line != '':
    headers.append(line)
    seqs.append(f.readline())
    print len(seqs[-1])
    line = f.readline()

seqs_new = [''] * len(seqs)
for i in range(len(seqs[0])):
    found_gap = False
    for j in range(len(seqs)):
        if seqs[j][i] == '-':
            found_gap = True
            break
    if not found_gap:
        for j in range(len(seqs)):
            seqs_new[j] += seqs[j][i]

for j in range(len(seqs)):
    fout.write(headers[j])
    fout.write(seqs_new[j])

f.close()
fout.close()
        
