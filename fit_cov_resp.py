#!/usr/bin/env python3
# A small script to fit the resp charge for covalent system
# by Junhao Li
import sys
import os
import math

global SPR, ACE, NME, AMBTYPE
#SPR: special residues, deprotonated with new charges
SPR = {'CYS':'CYI','GLU':'GUU','ASP':'APP','LYS':'LYY',
       'HIS':'HII','HIP':'HII','HIE':'HII','HID':'HII',
       'SER':'SRR'}
#ACE and NME (NMA in Schrodinger), amber+schrodinger names
#The shrodinger names are taken from neighboring residues
ACE = {'CH3':-0.3662,'H1':0.1123,'H2':0.1123,'H3':0.1123,
       'C':0.5972,'O':-0.5679,'CA':-0.3662,'HA':0.1123,
       'HA2':0.1123, 'HA3':0.1123, '1HH3':0.1123,
       '2HH3':0.1123, '3HH3':0.1123}
NME = {'N':-0.4157,'H':0.2719,'C':-0.1490,'CA':-0.1490,
       'H1':0.0976,'H2':0.0976,'H3':0.0976,'HA':0.0976,
       'HA2':0.0976,'HA3':0.0976, 'CH3':-0.1490,
       '1HH3':0.0976, '2HH3':0.0976, '3HH3':0.0976}
# this is just for CYS, more residue type can be import from amino12.lib
# and change AMBTYPE to a list of different residue dictionaries
AMBTYPE = {'N':'N','H':'H','C':'C','O':'O','CA':'CX','CB':'2C',
           'SG':'SH','HB2':'H1','HB3':'H1','HA':'H1'}

#def replace_cap_mc(line,atmNam):
#    if atmNam in ['N','C']:
#        new = atmNam+'1'
#        return line[0:12]+line[12:16].replace(atmNam+' ',new)+line[16:]
#    else:
#        return line

def repalce_mol2_type(inp,resid):
    ante = "%7d %-8s %10.4lf %10.4lf %10.4lf %-6s %5d %-8s %9.6lf\n"
    flag = False
    out = []
    for line in open(inp,'r'):
        if 'TRIPOS>ATOM' in line:   flag = True
        elif 'TRIPOS>BOND' in line: flag = False
        if flag and 'ATOM' not in line:
            atm = line.strip().split()
            if atm[6] == str(resid):
                newline = ante % (int(atm[0]),atm[1],float(atm[2]),float(atm[3]),
                float(atm[4]), AMBTYPE[atm[1]], int(atm[6]), atm[7],float(atm[8]))
                out.append(newline)
            else: out.append(line)
        else:
            out.append(line)
    return out

def get_init_chg(pdb):
    # return a list masking the positions of capped atoms
    # read the REMARK of cov bond manully added to the _opt.pdb by psi4_esp.py
    # and split the pdb into reactive residue(s) and ligand
    ind, out, ress, resi, resn = 0, [], {}, {}, []
    resi['LIG'] = []
    lig = open('ligonly.pdb','w')
    new = open('ligress.pdb','w')
#    cap = []
    n_cap = 0
    for line in open(pdb, 'r'):
        if 'REMARK' in line: bondinfo = line.strip().replace('REMARK', ' ')
        if line[0:6] in ['ATOM  ', 'HETATM']:
            atmNam = line[12:16].strip()
            resNam = line[17:20].strip()
            resNum = line[22:26].strip()
            if resNam != 'ACE' and resNam != 'NME':
                new.write(line[0:17].replace('HETATM','ATOM  ')+'LIS'+line[20:])
            if resNam == 'ACE':
                out.append(ACE[atmNam])
                n_cap += 1
#                cap.append(replace_cap_mc(line,atmNam))
            elif resNam in ['NME','NMA']:
                out.append(NME[atmNam])
#                cap.append(replace_cap_mc(line,atmNam))
                n_cap += 1
            else:
                out.append(0.0000)
            if resNam in SPR.keys():
                if SPR[resNam] not in ress.keys():
                    ress[SPR[resNam]] = []
                    resi[SPR[resNam]] = []
                    resn.append(resNum)
                ress[SPR[resNam]].append(line.replace(resNam,SPR[resNam]))
                resi[SPR[resNam]].append(ind)
            if resNam in ['LIG', 'MOL', 'UNK']:
                lig.write(line)
                resi['LIG'].append(ind)
            ind += 1
    new.close()
    lig.close()
    ires = 0
    # Note that the CYI is a residue now, its mol2 can be written
    for res in ress.keys():
        with open('res{}.pdb'.format(str(ires)),'w') as resF:
            resF.write(''.join(ress[res]))
#            resF.write(''.join(cap))
        ires += 1
    print("number of capped atoms are {}".format(n_cap))
    return [out,resi,ress.keys(),resn,bondinfo]

def write_init_chg(chgi, fout='qin.tmpl'):
    fp = open(fout, 'w')
    for i in range(len(chgi)):
        j = i + 1
        fp.write("{:10.6f}".format(float(chgi[i])))
        if float("{:.1f}".format(math.fmod(j,8))) == 0.0:
            fp.write('\n')
    fp.close()
    return 0

def mod_resp_inp(chgi,respin,stage=1,fout='resp1.in'):
    # Add a '-99' constraint to the capped atoms, because their
    # charges are taken from the amber force field directly
    i, j = -2, 0
    flag = False
    fp = open(fout,'w')
    for line in open(respin,'r'):
        j += 1
        if "Resp charges for" in line and j > 1:
            flag = True
            i += 1
            fp.write(line)
            continue
        if not flag:
            fp.write(line)
            if 'qwt' in line and stage == 1:
                fp.write(' iqopt = 2,\n')
        elif i == -2:
            fp.write(line)
        elif i == -1:
            fp.write(line)
            i += 1
            continue
        elif len(line.split()) > 1:
            #print(i)
            atmic = int(line.split()[0])
            state = line.split()[1]
            if chgi[i] == 0.0:
                fp.write(line)
            else:
                fp.write('{:5d}{:5d}\n'.format(atmic,-99))
            i += 1
        else:
            fp.write(line)
    fp.close()
    return 0

if len(sys.argv) < 3:
    print("usage: {} mol.pdb mol.esp".format(sys.argv[0]))
    print('put mol.pdb, mol.esp into current directory')
    exit(1)

pdb, esp = sys.argv[1], sys.argv[2]
#pdb, esp = 'mol.pdb', 'mol.esp'
base_pdb = '.'.join(pdb.split('.')[0:-1])
base_esp = '.'.join(esp.split('.')[0:-1])
# get_init_chg will also unify and split the original MAE exported pdb:
# ligress.pdb, ligonly.pdb, res0.pdb, res1.pdb, ... 
pdb_info = get_init_chg(pdb=pdb)
init_chg = pdb_info[0]
write_init_chg(init_chg,fout='cap.chg')
print(init_chg)
print(sum(init_chg),len(init_chg))

if not os.path.isdir('tmp'): os.mkdir('tmp')
os.system(f"$AMBERHOME/bin/espgen -i {esp} -o {base_esp}.dat")
os.chdir('tmp')
cmd0 = f"antechamber -fi gesp -i ../{esp} -fo ac -o {base_esp}.ac -c resp"
os.system(cmd0)
os.chdir('..')

mod_resp_inp(init_chg,'tmp/ANTECHAMBER_RESP1.IN',stage=1,fout='resp1.in')
mod_resp_inp(init_chg,'tmp/ANTECHAMBER_RESP2.IN',stage=2,fout='resp2.in')
cmd1 = f"$AMBERHOME/bin/resp -O -i resp1.in -o resp1.out -p resp1.pch -t resp1.chg -q cap.chg -e {base_esp}.dat"
cmd2 = f"$AMBERHOME/bin/resp -O -i resp2.in -o resp2.out -p resp2.pch -t resp2.chg -q resp1.chg -e {base_esp}.dat"
print(f"running custom resp fit\n{cmd1}\n{cmd2}\n")
os.system(cmd1)
os.system(cmd2)
# now the resp2.chg can be used, the sum charges of the capped atoms should zero.

#TODO fix the cases with more than one covalent bound residues
# 
deprot = [x for x in pdb_info[2] if x != 'LIG']
resi = pdb_info[1]
ress = {'ligress':'LIS','ligonly':'LIG','res0':deprot[0]}
all_resp = []
for line in open('resp2.chg','r'):
    all_resp += [float(x) for x in line.strip().split()]

if not os.path.isfile('ligress.chg'): os.link('resp2.chg', 'ligress.chg')
ligonly = [all_resp[x] for x in resi[ress['ligonly']]]
res0    = [all_resp[x] for x in resi[ress['res0']]]
print(f"ligand atoms: {len(ligonly)}\n{deprot[0]} atoms:  {len(res0)}\nTotal atoms:  {len(all_resp)}")
print(f"ligand charge {sum(ligonly)}\n{deprot[0]} charge: {sum(res0)}\nTotal charge: {sum(all_resp)}")
write_init_chg(ligonly,fout='ligonly.chg')
write_init_chg(res0, fout='res0.chg')

for prefix in ['ligress','ligonly','res0']:
    res = ress[prefix]
    atm, dr = 'gaff2', '  '
    if 'lig' not in prefix: atm='amber'
    if prefix == 'ligress': dr='-dr no'
    cmd3 = f"$AMBERHOME/bin/antechamber -fi pdb -i {prefix}.pdb -c rc -cf {prefix}.chg -fo ac -o {prefix}.ac -at {atm} {dr}"
    cmd4 = f"$AMBERHOME/bin/antechamber -fi pdb -i {prefix}.pdb -c rc -cf {prefix}.chg -fo mol2 -o {prefix}.mol2 -at {atm} {dr}"
    cmd5 = f"$AMBERHOME/bin/prepgen -i {prefix}.ac -o {prefix}.prep -rn {res}"
    cmd6 = f"$AMBERHOME/bin/parmchk2 -i {prefix}.mol2 -o {prefix}.frcmod -f mol2 -a Y"
    print(f"running:\n{cmd3}\n{cmd4}\n{cmd5}\n{cmd6}\n")
    os.system(cmd3)
    os.system(cmd4)
    if prefix == 'ligonly': os.system(cmd5)
    else:
        newlines = repalce_mol2_type(f'{prefix}.mol2',pdb_info[3][0])
        with open(f'{prefix}.mol2','w') as newmol2:
            newmol2.write(''.join(newlines))
    if prefix == 'ligress':
        with open('tmp.in','w') as fp:
            fp.write(f'p = loadmol2 {prefix}.mol2\n')
            fp.write(f'bond {pdb_info[4]}\n')
            fp.write(f'savemol2 p {prefix}.mol2 1\nquit\n')
        os.system('tleap -f tmp.in')
    os.system(cmd6)

FF='source leaprc.protein.ff14SB'
res0=deprot[0]
with open('res0.in','w') as leap1:
    leap1.write(f'''{FF}
{res0} = loadmol2 res0.mol2
list
desc {res0}
set {res0} head {res0}.{res0}.N
set {res0} tail {res0}.{res0}.C
set {res0} restype protein
desc {res0}
saveoff {res0} res0.lib
quit''')

os.system("$AMBERHOME/bin/tleap -f res0.in")
# for better visualization:
os.system('antechamber -fi mol2 -i ligress.mol2 -fo mol2 -o ligress_sybyl.mol2 -at sybyl -dr no')
print("Files to be used for covMD: res0.lib, ligonly.prep ligress.frcmod")
print("Files to be used for deprotonated CYS MD: res0.lib, ligonly.prep, ligonly.frcmod")
