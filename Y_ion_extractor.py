usage = \
"""
Usage: python Y_ion_extractor D:/xx/xx/pGlyco3.cfg
Output: extract_Y_ions.txt in 'output_dir' pointed in pGlyco3.cfg
"""

import os
import struct
import pandas as pd
import numpy as np

mass_proton = 1.007276
mass_isotope = 1.0033
tol = 20.0
ppm = True
deisotope=False

glyco_mass = {}
glyco_mass['N'] = 203.079372533
glyco_mass['H'] = 162.0528234315
glyco_mass['F'] = 146.0579088094
glyco_mass['A'] = 291.09541652769997
glyco_mass['G'] = 307.09033114979997
glyco_mass['X'] = 132.0422587452
glyco_mass['pH'] = 242.0191544315
glyco_mass['aH'] = 179.0793635


Yions = ['N(1)F(1)', 'N(2)F(1)']
pGlycoResult = "pGlycoDB-GP-FDR-Pro-Quant-Site.txt"

def calc_glycan_mass(glycan):
    glycos = glycan.strip(')').split(')')
    mass = 0
    for glyco in glycos:
        glyco, n = glyco.split('(')
        mass += glyco_mass[glyco]*int(n)
    return mass
    
Yion_masses = np.array([calc_glycan_mass(glycan) for glycan in Yions])

class pf2reader:
    def __init__(self, pf2=None):
        self.scan2spec = {}
        self.pf2 = None
        if pf2:
            self.open(pf2)
        
    def close(self):
        if self.pf2 is not None:
            self.pf2.close()
            self.pf2 = None
        
    def open(self, pf2):
        self.close()
        self.pf2 = open(pf2,'rb')
        
        self.scan2spec = {}
        with open(pf2+'idx','rb') as f:
            while True:
                chunk = f.read(8)
                if not chunk: break
                scan, index = struct.unpack('2i',chunk)
                self.scan2spec[scan] = index
                
    def read_peaklist(self, scan):
        # only read the (mz, inten) list
        self.pf2.seek(self.scan2spec[scan])
        _scan, nPeak = struct.unpack("2i",self.pf2.read(8))
        mz_int = struct.unpack(str(nPeak*2)+"d", self.pf2.read(nPeak*2*8))
        masses = []
        intens = []
        for i in range(nPeak):
            mz = mz_int[i*2]
            inten = mz_int[i*2+1]
            masses.append(mz)
            intens.append(inten)
        return np.array(masses), np.array(intens)

def read_until(file, until):
    lines = []
    while True:
        line = file.readline().strip()
        if line.startswith(until):
            break
        else:
            lines.append(line)
    return lines

def find_line(lines, start):
    for line in lines:
        if line.startswith(start):
            return line
    return None

def parse_scan_from_TITLE(pfind_title):
    return int(pfind_title.split('.')[-4])

class MGFReader():
    def __init__(self, mgf=None):
        self.scan2spec = {}
        if mgf is not None:
            self.open(mgf)

    def close(self):
        pass

    def open(self, mgf):
        self.scan2spec = {}
        print("Loading MGF ...")
        with open(mgf) as f:
            while True:
                line = f.readline()
                if not line: break
                if line.startswith('BEGIN IONS'):
                    lines = read_until(f, 'END IONS')
                    masses = []
                    intens = []
                    scan = None
                    for line in lines:
                        if line[0].isdigit():
                            mass,inten = [float(i) for i in line.strip().split()]
                            masses.append(mass)
                            intens.append(inten)
                        elif line.startswith('SCAN='):
                            scan = int(line.split('=')[1])
                    if not scan:
                        title = find_line(lines, 'TITLE=')
                        scan = parse_scan_from_TITLE(title)
                    self.scan2spec[scan] = (np.array(masses), np.array(intens))
                    print(f'Scan={scan}', end='\r')
        print("Finish loading MGF")

    def read_peaklist(self, scan):
        return self.scan2spec[scan]

def get_ms2_reader(pf2):
    if os.path.isfile(pf2):
        return pf2reader(pf2)
    else:
        return MGFReader(pf2[:-3]+'mgf')

def match(masses, intens, query_masses):
    query_masses = np.sort(query_masses)
    if ppm: tols = query_masses*tol*1e-6
    else: tols = np.ones_like(query_masses)*tol
    i = 0
    j = 0
    query_intens = np.zeros_like(query_masses)
    while i < len(masses) and j < len(query_masses):
        if masses[i] > (query_masses[j]+tols[j]):
            j += 1
        elif masses[i] < (query_masses[j]-tols[j]):
            i += 1
        else:
            query_intens[j] += intens[i]
            i += 1
    return query_intens
    
if __name__ == "__main__":
    import sys
    if len(sys.argv) == 1:
        print(usage)
        sys.exit(-1)
    cfg = sys.argv[1]
    
    with open(cfg) as f:
        lines = f.readlines()
        raw_names = []
        pf2_readers = {}
        for line in lines:
            if line.startswith('output_dir'):
                output_dir = line[line.find('=')+1:].strip()
            elif line.startswith('file'):
                _dir, filename = os.path.split(line[line.find('=')+1:].strip())
                if filename.lower().endswith('.raw'):
                    raw_name = filename[:-4]
                else:
                    raw_name = filename[:filename.rfind('_')]
                raw_names.append(raw_name)
                pf2_readers[raw_name] = get_ms2_reader(os.path.join(_dir, raw_name+"_HCDFT.pf2"))
                
                
    df = pd.read_csv(os.path.join(output_dir, pGlycoResult), sep='\t')
    
    for _head in df.columns.values:
        if _head.startswith('Glycan('):
            glycan_head = _head
            glycan_list = glycan_head[len('Glycan('):-1].split(',')
            break
            
    def glycan_composition_to_vector(glycan_comp, glycan_list):
        ret = np.zeros(len(glycan_list))
        glycos = glycan_comp.strip(')').split(')')
        for glyco in glycos:
            glyco, n = glyco.split('(')
            if glyco in glycan_list:
                idx = glycan_list.index(glyco)
                ret[idx] = float(n)
        return ret
        
    Yion_vectors = [glycan_composition_to_vector(glycan_comp, glycan_list) for glycan_comp in Yions]
    
    matched_intens = []
    for rawname, scan, charge, peptide_mass, glycan_vector in df[f'RawName;Scan;Charge;PeptideMH;{glycan_head}'.split(';')].values:
        print(rawname, scan)
        peptide_mass -= mass_proton
        glycan_vector = np.array([float(g) for g in glycan_vector.strip().split(' ')])
        masses, intens = pf2_readers[rawname].read_peaklist(scan)
        
        extract_intens = []
        for Ymass,Yvector in zip(Yion_masses,Yion_vectors):
            if np.all(glycan_vector>=Yvector):
                
                ions = (Ymass + peptide_mass)/np.arange(1, charge) + mass_proton
                ion_intens = match(masses, intens, ions)
                
                if deisotope:
                    preisotope = ions - mass_isotope/np.arange(1, charge)
                    preisotope = match(masses, intens, preisotope)
                    
                    ion_intens = ion_intens[preisotope==0]
                
                extract_intens.append(np.sum(ion_intens))
            else:
                extract_intens.append(0)
            
        matched_intens.append(extract_intens)
    matched_intens = np.array(matched_intens)
    
    df[['Y-'+ion for ion in Yions]] = matched_intens
    
    df.to_csv(os.path.join(output_dir, 'extract_Y_ions.txt'), sep='\t', index=False)
    print("result saved as '%s'"%os.path.join(output_dir, 'extract_Y_ions.txt'))
