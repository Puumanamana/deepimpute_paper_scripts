import sys
from os.path import splitext
import numpy as np

filename = sys.argv[1]
_,method,ncells,trial = splitext(filename)[0].split('_')

n_lines = sum(1 for _ in open(filename))

handler = open(filename)

memory = []
for line in handler:
    data = line.strip().split(" ")
    data = [ x for x in data if len(x)>0]

    try:
        mem = int(data[3])
    except:
        continue
    memory.append(mem)
handler.close()

if len(memory) == 0:
    print("Empty memory vector. Aborting")
    exit(1)

memory = np.array(memory)
peak_usage = memory[0] - memory.min()

print("{} - Check: {} memory values ({} lines)".format(filename,len(memory),n_lines))
print("{}: Peak memory usage = {} MB".format(filename,peak_usage/1e3))

with open('summary.csv','a') as handle:
    handle.write(','.join([method,ncells,trial,str(peak_usage),"peakMem"]))
    handle.write('\n')
    
