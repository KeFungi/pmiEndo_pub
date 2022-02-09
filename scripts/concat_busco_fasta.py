#!/usr/bin/python
import os
import sys

folder = sys.argv[1] + '/'
splist_file = sys.argv[2]

filelist = [folder + f for f in os.listdir(folder) if f.endswith("_ali.fasta")]
splist = []
with open(splist_file, 'r') as spfile:
    for line in spfile:
        splist.append(line.strip())


data = dict.fromkeys(splist, '')


for f in filelist:
    for line in open(f, 'r'):
        if line.startswith('>'):
            key=line.strip()[1:]
            if key in data:
                pass
            else:
                print("error: unknown species code in fasta file")
        else:
            data[key] += line.strip()
# output
with open('fungi_supermatrix_ali.fasta', 'w') as out:
    for key in data:
        out.write('>' + key + '\n')
        out.write(data[key] + '\n')

