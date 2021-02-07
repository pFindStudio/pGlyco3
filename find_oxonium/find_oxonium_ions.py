scan_txt = r"C:\workspace\pf2_alphapept-branch\find_oxonium_ions\scan.txt"
pf2 = r"C:\DataSets\test_AP\Thermo\hela500ng\HeLa_500ng_HCDFT.pf2"
output = r"C:\workspace\pf2_alphapept-branch\find_oxonium_ions\oxonium.txt"
max_mz = 300
min_rel_inten = 0.1
bin_per_Da = 0.1

import struct
import numpy as np

scan_set = set()
with open(scan_txt) as f: scan_set = set([int(scan.strip()) for scan in f.readlines()])

scanidx = {}
with open(pf2+'idx','rb') as f:
    while True:
        chunk = f.read(8)
        if not chunk: break
        scan, index = struct.unpack('2i',chunk)
        scanidx[scan] = index

pf2 = open(pf2, 'rb')

def read_a_peaklist(scan):
    # only read the (mz, inten) list
    pf2.seek(scanidx[scan])
    scan, nPeak = struct.unpack("2i",pf2.read(8))
    mz_int = struct.unpack(str(nPeak*2)+"d", pf2.read(nPeak*2*8))
    masses = []
    intens = []
    for i in range(nPeak):
        mz = mz_int[i*2]
        inten = mz_int[i*2+1]
        masses.append(mz)
        intens.append(inten)
    return np.array(masses), np.array(intens)
        
def find_oxonium_ions(scan):
    masses, intens = read_a_peaklist(scan)
    max_inten = max(intens)
    bins = np.zeros(int(max_mz/bin_per_Da))
    for mz, inten in zip(masses, intens):
        idx = int(mz/bin_per_Da)
        if idx >= len(bins): break 
        if inten/max_inten >= min_rel_inten:
            bins[idx] = 1
    return bins
    
counts = np.zeros(int(max_mz/bin_per_Da))
for scan in list(scan_set):
    if scan not in scanidx:
        scan_set.remove(scan)
        continue
    print(scan)
    counts += find_oxonium_ions(scan)
frequency = counts/len(scan_set)

with open(output,'w') as f:
    f.write('mz\tfrequency\n')
    for mz,freq in enumerate(frequency):
        f.write('{}\t{}\n'.format(mz*bin_per_Da, freq))

pf2.close()
