from chemical import *
import sys

if __name__ == "__main__":
    # marker 138 formula: C(7)H(7)N(1)O(2)
    if len(sys.argv) == 1:
        print("Usage: python [-f formula] [-ini inifile -formula_pos n -mass_pos m] [-r replacement]")
        print("marker 138 formula: python -f C(7)H(7)N(1)O(2) -r N~15N")
        sys.exit(-1)
    argd = {}
    for i in range(1, len(sys.argv), 2):
        argd[sys.argv[i]] = sys.argv[i+1]
    if '-r' in argd:
        replacement = [item.split('~') for item in argd['-r'].split(',')]
    else:
        replacement = []
    if '-f' in argd:
        formula = argd['-f']
        if replacement:
            f, m = replace_element_and_calc_mass(formula, replacement)
            print(f, "%.10f"%m)
        else:
            print(str(calc_formula_mass(formula)))
    elif '-ini' in argd:
        formula_pos = int(argd['-formula_pos'])
        mass_pos = int(argd['-mass_pos'])
        with open(argd['-ini']) as f: lines = f.readlines()
        with open(argd['-ini']+".bak", "w") as f:
            for line in lines:
                print(line)
                items = line.split(" ")
                if replacement:
                    formula, mass = replace_element_and_calc_mass(items[formula_pos], replacement)
                    items[formula_pos], item[mass_pos] = formula, str(mass)
                else:
                    mass = calc_formula_mass(items[formula_pos].strip())
                    items[mass_pos] = str(mass)
                f.write(" ".join(items))
