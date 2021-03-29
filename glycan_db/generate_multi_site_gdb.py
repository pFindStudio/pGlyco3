import sys
import itertools

max_comb=2
max_len=25

gdb = sys.argv[1]

def canon_to_comp(canon, glyco_list):
    comp = [0]*len(glyco_list)
    items = canon.strip('(').split('(')
    for item in items:
        item = item.strip(')')
        idx = glyco_list.index(item)
        comp[idx] += 1
    return tuple(comp)

with open(gdb) as f:
    head = f.readline().strip()
    glycans = [line.strip() for line in f.readlines()]

glyco_list = head.split(',')
ret = [glycan for glycan in glycans]
comp_set = set([canon_to_comp(glycan, glyco_list) for glycan in ret])
glycans.append("")
glycans.sort(key=lambda x: (len(x), x))

with open(gdb[:-4]+'-multi.gdb','w') as f:
    f.write(head + '\n')
    for comb in itertools.combinations_with_replacement(glycans, max_comb):
        new_glycan = ''.join(comb)
        if not new_glycan: continue
        new_comp = canon_to_comp(new_glycan, glyco_list)
        if new_comp not in comp_set and sum(new_comp) <= max_len:
            ret.append(new_glycan)
            comp_set.add(new_comp)
            print(f'{len(ret)} structures generated')
    ret.sort(key=lambda x: (len(x), x))
    f.write('\n'.join(ret)+'\n')
