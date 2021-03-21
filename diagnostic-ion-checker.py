import os
import struct
import numba
import pandas as pd
import numpy as np

def match(masses, intens, marker, tol):
    i = 0
    inten = 0
    while i < len(masses):
        if masses[i] > (marker+tol):
            break
        elif masses[i] < (marker-tol):
            i += 1
        else:
            inten += intens[i]
            i += 1
    return inten

def check_marker_pf2(pf2_list, marker_list):
    tol = 0.02

    raw_list = []
    scan_list = []
    matched_list = []
    for pf2 in pf2_list:
        pf2idx = {}
        with open(pf2+'idx','rb') as f:
            while True:
                chunk = f.read(8)
                if not chunk: break
                scan_no, index = struct.unpack('2i',chunk)
                pf2idx[scan_no] = index
        
        raw_name = os.path.basename(pf2[:pf2.rfind('_')]).upper()
        print(raw_name)

        f = open(pf2, 'rb')
        for scan, idx in pf2idx.items():
            f.seek(idx)
            scan, nPeak = struct.unpack("2i",f.read(8))
            mz_int = np.array(struct.unpack("%dd"%(nPeak*2), f.read(nPeak*2*8)))
            masses = mz_int[0::2]
            intens = mz_int[1::2]
            max_inten = np.max(intens)
            matched_inten = []
            for marker in marker_list:
                inten = match(masses, intens, marker, tol)
                matched_inten.append(inten/max_inten)
            raw_list.append(raw_name)
            scan_list.append(scan)
            matched_list.append(matched_inten)
        f.close()
    
    matched_list = np.array(matched_list)
    df = pd.DataFrame(None)
    df['RawName'] = raw_list
    df['Scan'] = scan_list
    for i,marker in enumerate(marker_list):
        df[f'Marker_{marker:.2f}'] = matched_list[:,i]
    return df
    
if __name__ == "__main__":
    argd = {}
    import sys
    for i in range(1, len(sys.argv), 2):
        argd[sys.argv[i]] = sys.argv[i+1]
    
    pf2_list = []
    if os.path.isdir(argd['-pf2']):
        for pf2 in os.listdir(argd['-pf2']):
            if pf2.endswith('.pf2'):
                pf2_list.append(os.path.join(argd['-pf2'], pf2))
        outdir = argd['-pf2']
    else:
        pf2_list.append(argd['-pf2'])
        outdir = os.path.split(argd['-pf2'])[0]
        
    marker_list = [eval(i) for i in argd['-marker'].split(';')]
    
    df = check_marker_pf2(pf2_list, marker_list)
    df.to_csv(os.path.join(outdir, 'marker-ion.csv'))
