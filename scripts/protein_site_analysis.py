import os

glycan_fdr_thres = 0.01
peptide_fdr_thres = 0.1

def ProteinAnalysis(txt):
    with open(txt) as f:
        head = f.readline().strip().split('\t')
        headidx = dict(zip(head, range(len(head))))
        def get(items, name):
            return items[headidx[name]]
            
        glycan_head = "Glycan"
        for name in head:
            if name.startswith('Glycan('):
                glycan_head = name
                break
        
        protein_sites = []
        lines = f.readlines()
        for line in lines:
            items = line.strip().split('\t')
            glycan_fdr = get(items, 'GlycanFDR')
            peptide_fdr = get(items, 'PeptideFDR')
            total_fdr = get(items, 'TotalFDR')
            if float(glycan_fdr) > glycan_fdr_thres: continue
            if float(peptide_fdr) > peptide_fdr_thres: continue
            spec = get(items, 'GlySpec')
            raw = get(items, 'RawName')
            scan = get(items, 'Scan')
            ETDscan = get(items, 'ETDScan')
            peptide = get(items, 'Peptide')
            pepsite = int(get(items, 'GlySite'))
            mod = get(items, 'Mod')
            glycan = get(items, glycan_head)
            mono = get(items, 'MonoArea')
            isotope = get(items, 'IsotopeArea')
            proteins = get(items, 'Proteins').split(';')
            genes = get(items, 'Genes').split(';')
            prosites = get(items, 'ProSites').split(';')
            SiteGroups = get(items, 'LocalizedSiteGroups')
            if not SiteGroups:
                for i,protein in enumerate(proteins):
                    gene = genes[i]
                    prositegroup = peptide[pepsite-1] + prosites[i] + ':' + peptide[pepsite-1] + prosites[i]
                    protein_sites.append([protein, gene, prositegroup, prosites[i], prosites[i], "", "-1", "", "", "", mono, isotope, peptide, glycan, prosites[i], mod, spec, raw, scan, ETDscan, glycan_fdr, peptide_fdr, total_fdr])
            else:
                SiteGroups = [group.split(',') for group in SiteGroups.split(';')]
                for s1,s2,site_glycan,prob in SiteGroups:
                    site1 = int(s1[1:])
                    site2 = int(s2[1:])
                    site_glycan = site_glycan[1:-1]
                    prob = prob
                    for i,protein in enumerate(proteins):
                        gene = genes[i]
                        prosite = int(prosites[i])
                        prosite1 = str(prosite+site1-pepsite)
                        prosite2 = str(prosite+site2-pepsite)
                        prositegroup = s1[0] + prosite1 + ':' + s2[0] + prosite2
                        protein_sites.append([protein, gene, prositegroup, prosite1, prosite2, site_glycan, prob, "1" if prosite1==prosite2 else "0", "1" if prosite1!=prosite2 else "0", s1+':'+s2, mono, isotope, peptide, glycan, str(prosite), mod, spec, raw, scan, ETDscan, glycan_fdr, peptide_fdr, total_fdr])
    protein_sites.sort(key = lambda x: (x[0], int(x[3]), int(x[4]), x[5], -float(x[6])))
    protein_sites.insert(0, ["Protein", "Gene", "ProteinSiteGroup", "StartSite", "EndSite", 'Localized'+glycan_head, "SiteProbability", "IsUniqueSite", "IsGroupSite", "PeptideSiteGroup", "MonoArea", "IsotopeArea", "Peptide", "AllSiteGlycan", "ProSite", "Mod", "Spectrum", "RawName", "Scan", "ETDScan", 'GlycanFDR', 'PeptideFDR', 'TotalFDR'])
    return protein_sites

def ProteinOutput(txt, protein_sites):
    with open(os.path.join(os.path.dirname(txt), 'pGlyco3-ProSite-Analysis.txt'),'w') as f:
        f.writelines(['\t'.join(items)+'\n' for items in protein_sites])
    print('results saved as "%s"'%os.path.join(os.path.dirname(txt), 'pGlyco3-ProSite-Analysis.txt'))
        
if __name__ == '__main__':
    import sys
    if len(sys.argv) != 2:
        print("Usage: python protein_site_analysis.py D:/xx/xx/pGlycoDB-GP-FDR-Pro-Quant-Site.txt")
        sys.exit(-1)
    protein_sites = ProteinAnalysis(sys.argv[1])
    ProteinOutput(sys.argv[1], protein_sites)
                    
