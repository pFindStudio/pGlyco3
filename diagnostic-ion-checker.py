usage = \
"""
Usage:
    python diagnostic-ion-checker.py -mgf mgf_containing_folder -marker 204.087,138.055
Output:
    marker-ion.csv in the mgf_containing_folder

Some candidate marker ions:
    | Marker  | Glyco            | Formula w/o proton |
    | --------| ---------------- | ------------------ |
    | 109.028 | Hex fragment     | C6H4O2             |
    | 115.039 | Hex fragment     | C5H6O3             |
    | 126.055 | HexNAc fragment  | C6H7O2N1           |
    | 127.039 | Hex–H2O*2        | C6H6O3             |
    | 138.055 | HexNAc fragment  | C7H7O2N1           |
    | 144.066 | HexNAc fragment  | C6H9O3N1           |
    | 163.060 | Hex              | C6H10O5            |
    | 168.066 | HexNAc–H2O*2     | C8H9O3N1           |
    | 186.076 | HexNAc–H2O       | C8H11O4N1          |
    | 204.087 | HexNAc           | C8H13O5N1          |
    | 274.092 | NeuAc–H2O        | C11H15O7N1         |
    | 290.087 | NeuGc–H2O        | C11H15O8N1         |
    | 292.103 | NeuAc            | C11H17O8N1         |
    | 308.098 | NeuGc            | C11H17O9N1         |
    | 366.140 | Hex+HexNAc       | C14H23O10N1        |
    | 657.140 | Hex+HexNAc+NeuAc | C25H40O18N2        |
    | 673.230 | Hex+HexNAc+NeuGc | C25H40O19N2        |
"""

import os
import struct
import numba
import pandas as pd
import numpy as np

from Y_ion_extractor import get_ms2_reader

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
    basepeak_list = []
    for pf2 in pf2_list:
        raw_name = os.path.basename(pf2[:pf2.rfind('_')]).upper()
        print(raw_name)
        ms_reader = get_ms2_reader(pf2)

        for scan in ms_reader.scan2spec.keys():
            masses, intens = ms_reader.read_peaklist(scan)
            matched_inten = []
            for marker in marker_list:
                inten = match(masses, intens, marker, tol)
                matched_inten.append(inten)
            raw_list.append(raw_name)
            scan_list.append(scan)
            matched_list.append(matched_inten)
            basepeak_list.append(np.max(intens))
    
    matched_list = np.array(matched_list)
    df = pd.DataFrame(None)
    df['RawName'] = raw_list
    df['Scan'] = scan_list
    df['BasePeak'] = basepeak_list
    for i,marker in enumerate(marker_list):
        df[f'Marker_{marker:.2f}'] = matched_list[:,i]
    return df
    
if __name__ == "__main__":
    argd = {}
    import sys
    if len(sys.argv) == 1:
        print(usage)
        sys.exit(-1)
    for i in range(1, len(sys.argv), 2):
        argd[sys.argv[i]] = sys.argv[i+1]
    
    pf2_list = []
    if os.path.isdir(argd['-mgf']):
        for mgf in os.listdir(argd['-mgf']):
            if mgf.endswith('.mgf'):
                pf2_list.append(os.path.join(argd['-mgf'], mgf[:-3]+'pf2'))
            elif mgf.endswith('.pf2'):
                pf2_list.append(os.path.join(argd['-mgf'], mgf))
        outdir = argd['-mgf']
    else:
        if argd['-mgf'].endswith('.mgf'):
            pf2_list.append(argd['-mgf'][:-3]+'pf2')
        elif argd['-mgf'].endswith('.pf2'):
            pf2_list.append(argd['-mgf'])
        outdir = os.path.split(argd['-mgf'])[0]
        
    marker_list = [eval(i) for i in argd['-marker'].split(',')]
    
    df = check_marker_pf2(pf2_list, marker_list)
    df.to_csv(os.path.join(outdir, 'marker-ion.csv'))
