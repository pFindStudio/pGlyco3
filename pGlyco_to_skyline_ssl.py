import sys
import os

cfg = r"C:\folder\to\pglyco\result\pGlyco3.cfg"

raw_dict = {}
mod_dict = {}
mod_dict['Oxidation[M]'] = ('Oxidation (M)','15.994915')
mod_dict['Carbamidomethyl[C]'] = ('Carbamidomethyl (C)','57.021464')
mod_dict['Acetyl[AnyN-term]'] = ('Acetyl (N-term)','42.010565')
mod_dict['Acetyl[ProteinN-term]'] = ('Acetyl (N-term)','42.010565')

if len(sys.argv) == 2:
    cfg = sys.argv[1]
    lines = open(cfg).readlines()
    for line in lines:
        if line.startswith('file'):
            raw = line.split('=')[1].strip()
            if raw.lower().endswith('mgf'):
                raw = raw[:raw.rfind('_')]
            else:
                raw = raw[:-4]
            raw_dict[os.path.split(raw)[1]] = raw.replace('\\', '/')+'.raw'
        elif line.startswith('output_dir'):
            output_dir = line.split('=')[1].strip()
        elif line.startswith('pGlyco_type'):
            pGlycoType = line.split('=')[1].strip()
    pglyco = os.path.join(output_dir, '%s-GP-FDR.txt'%pGlycoType)
elif len(sys.argv) == 3:
    pglyco = sys.argv[1]
    raw_folder = sys.argv[2]
    output_dir = os.path.split(pglyco)[0].replace('\\','/')
else:
    print(
"""
Usage:
python pGlyco_to_skyline_ssl.py C:/output/folder/pGlyco.cfg
or 
python pGlyco_to_skyline_ssl.py C:/output/folder/pGlycoDB-GP-FDR.txt raw_folder
"""
    )
    sys.exit(-1)

skyline = os.path.join(output_dir, 'skyline')
try:
    os.makedirs(skyline)
except:
    pass

ssl = os.path.join(skyline,'pGlyco-skyline.ssl')
xml = os.path.join(skyline,'pGlyco-skyline.sky')
print(ssl)
with open(pglyco) as f, open(ssl,'w') as out:
    xml_start = """<?xml version="1.0" encoding="utf-8"?>
<srm_settings format_version="20.2" software_version="Skyline (64-bit) 20.2.0.286 (59e53e849)">
  <settings_summary name="Default">
    <peptide_settings>
      <enzyme name="Trypsin/P" cut="KR" no_cut="" sense="C" />
      <digest_settings max_missed_cleavages="3" />
      <peptide_prediction use_measured_rts="true" measured_rt_window="2" />
      <peptide_filter start="25" min_length="8" max_length="25" auto_select="true">
        <peptide_exclusions />
      </peptide_filter>
      <peptide_libraries pick="library" />
      <peptide_modifications max_variable_mods="3" max_neutral_losses="2">
        <static_modifications>
          <static_modification name="Carbamidomethyl (C)" aminoacid="C" formula="H3C2NO" unimod_id="4" short_name="CAM" />
          <static_modification name="Oxidation (M)" aminoacid="M" variable="true" formula="O" unimod_id="35" short_name="Oxi">
            <potential_loss formula="H4COS" massdiff_monoisotopic="63.998285" massdiff_average="64.10701" />
          </static_modification>"""
    xml_end = """
        </static_modifications>
        <heavy_modifications />
      </peptide_modifications></peptide_settings>
    <transition_settings>
      <transition_prediction precursor_mass_type="Monoisotopic" fragment_mass_type="Monoisotopic" optimize_by="None" />
      <transition_filter precursor_charges="2,3,4,5,6" product_charges="1" precursor_adducts="[M+H]" product_adducts="[M+]" fragment_types="p" small_molecule_fragment_types="f" fragment_range_first="m/z &gt; precursor" fragment_range_last="3 ions" precursor_mz_window="0" auto_select="true">
        <measured_ion name="N-terminal to Proline" cut="P" sense="N" min_length="3" />
      </transition_filter>
      <transition_libraries ion_match_tolerance="0.5" min_ion_count="0" ion_count="3" pick_from="all" />
      <transition_integration integrate_all="true" />
      <transition_instrument min_mz="50" max_mz="1500" mz_match_tolerance="0.055" />
      <transition_full_scan precursor_isotopes="Count" precursor_isotope_filter="5" precursor_mass_analyzer="centroided" precursor_res="10" retention_time_filter_type="ms2_ids" retention_time_filter_length="5">
        <isotope_enrichments name="Default">
          <atom_percent_enrichment symbol="H'">0.98</atom_percent_enrichment>
          <atom_percent_enrichment symbol="C'">0.995</atom_percent_enrichment>
          <atom_percent_enrichment symbol="C&quot;">0.99</atom_percent_enrichment>
          <atom_percent_enrichment symbol="N'">0.995</atom_percent_enrichment>
          <atom_percent_enrichment symbol="O&quot;">0.99</atom_percent_enrichment>
          <atom_percent_enrichment symbol="O'">0.99</atom_percent_enrichment>
          <atom_percent_enrichment symbol="Cl'">0.99</atom_percent_enrichment>
          <atom_percent_enrichment symbol="Br'">0.99</atom_percent_enrichment>
          <atom_percent_enrichment symbol="P'">0.99</atom_percent_enrichment>
          <atom_percent_enrichment symbol="S&quot;">0.99</atom_percent_enrichment>
          <atom_percent_enrichment symbol="S'">0.99</atom_percent_enrichment>
          <atom_percent_enrichment symbol="H&quot;">0.99</atom_percent_enrichment>
        </isotope_enrichments>
      </transition_full_scan>
    </transition_settings>
    <data_settings document_guid="e42b036a-613b-47cd-bf48-0bca60461e28" audit_logging="true" />
  </settings_summary>
</srm_settings>
"""
    xml_str = ""
    xml_set = set()
    def AddXml(glycan, glymass, aa, xml_set):
        if aa+glymass in xml_set: return ""
        xml_set.add(aa+glymass)
        mass = float(glymass)
        one_HexNAc = mass - 203.07937
        # return """
          # <static_modification name="{glymass}" aminoacid="{aa}" variable="true" massdiff_monoisotopic="{glymass}" massdiff_average="{glymass}" />""".format(glycan=glycan, glymass=glymass, aa=aa,one_HexNAc=one_HexNAc)
        return """
          <static_modification name="{glymass}" aminoacid="{aa}" variable="true" massdiff_monoisotopic="{glymass}" massdiff_average="{glymass}">
            <potential_loss massdiff_monoisotopic="{glymass}" massdiff_average="{glymass}" inclusion="Always" />
            <potential_loss massdiff_monoisotopic="{one_HexNAc}" massdiff_average="{one_HexNAc}" />
          </static_modification>""".format(glycan=glycan, glymass=glymass, aa=aa,one_HexNAc=one_HexNAc) 
        
    out.write('file	scan	charge	sequence	retention-time	glycan\n')
    head = f.readline().strip().split('\t')
    for i in head:
        if i.startswith('Glycan('):
            glycan_head = i
            break
    
    glynames = glycan_head[glycan_head.find('(')+1:-1].split(',')
    def get_glycan(glycan_str, glynames):
        glycos = [int(i) for i in glycan_str.strip().split(' ')]
        ret = ""
        for glyco, i in zip(glynames, glycos):
            if i > 0: ret += "%s(%d)"%(glyco,i)
        return ret
    
    headidx = dict(zip(head, range(len(head))))
    
    lines = f.readlines()
    for line in lines:
        items = line.strip().split('\t')
        raw = items[headidx['RawName']]
        scan = items[headidx['Scan']]
        rt = float(items[headidx['RT']])/60
        charge = items[headidx['Charge']]
        seq = items[headidx['Peptide']]
        seq = seq.replace('J','N')
        mod = items[headidx['Mod']]
        glymass = items[headidx['GlyMass']]
        glycan = items[headidx[glycan_head]]
        glycan = get_glycan(glycan, glynames)
        glysite = int(items[headidx['GlySite']])
        mod_list = [(glysite, (glycan, glymass))]
        
        aa = seq[glysite-1]
        xml_str += AddXml(glycan, glymass, aa, xml_set)
        
        if mod:
            mods = mod.strip(';').split(';')
            for mod in mods:
                site, mod = mod.split(',')
                mod_list.append((int(site), mod_dict[mod]))
            mod_list.sort(reverse=True)
        
        for site, (mod,mass) in mod_list:
            if site == 0: site = 1
            elif site == len(seq)+1: site = len(seq)
            else: site = site
            mod = '[%s]'%mass
            seq = seq[:site] + mod + seq[site:]
        if raw not in raw_dict:
            raw = os.path.join(raw_folder, raw+'.raw')
        else:
            raw = raw_dict[raw]
        out.write("{}\t{}\t{}\t{}\t{}\t{}\n".format(raw, scan, charge, seq, rt, glycan))
    with open(xml,'w') as xml:
        xml.write(xml_start+xml_str+xml_end)
    