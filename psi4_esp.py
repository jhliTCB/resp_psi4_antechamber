#!/public/home/lijunhao/soft/conda_envs/p4env/bin/python
#https://github.com/cdsgroup/resp/tree/master/examples
# resp package need to be patched
# by Junhao Li
import psi4
import resp
#from rdkit import Chem
#from rdkit.Chem import AllChem
import argparse

#def smi2xyz(mol):
#    mol = Chem.AddHs(mol)
#    AllChem.EmbedMolecule(mol, useExpTorsionAnglePrefs=True,useBasicKnowledge=True)
#    AllChem.UFFOptimizeMolecule(mol)
#    atoms = mol.GetAtoms()
#    string = "\n"
#    for i, atom in enumerate(atoms):
#        pos = mol.GetConformer().GetAtomPosition(atom.GetIdx())
#        string += "{} {} {} {}\n".format(atom.GetSymbol(), pos.x, pos.y, pos.z)
#    string += "units angstrom\n"
#    return string, mol

def read_par():
    parser = argparse.ArgumentParser(description="Calculation of electrophilicity index")
    parser.add_argument('-ismi',     dest='ismi',           type=str,  default='None',help="using a smile string as input")
    parser.add_argument('-ixyz',     dest='ixyz',           type=str,  default='None',help="using a xyz file as input")
    parser.add_argument('-ipdb',     dest='ipdb',        type=str,  default='None',  help='use a maestro/pymol exported PDB file as input')
    parser.add_argument('-no_opt',   dest='no_opt',        action='store_true',default=False,
        help='if -no_opt is given, then we will not perform optimization at B3LYP/6-31G*')
    parser.add_argument('-gau_temp', dest='gau',           action='store_true',default=False,help='generate a gaussian input template')
    parser.add_argument('-cm',       dest='chgmul',        type=str,  required=True,
        help='specify the charge and multiplicity, Required! Separated by comma, e.g. "-cm 0,1"; for negative charge, e.g. "-cm e1,1"')
    parser.add_argument('-c',        dest='core',          type=int,  default=8,help='specify the number of cores to use, default=8')
    parser.add_argument('-m',        dest='memory',        type=int,  default=16,help='''specify the number of memory (unit GB) to use, default=16;
           NOTE THAT, some basis sets will require huge amount of memory for property calculation!''')
    parser.add_argument('-dry',      dest='dry',       action='store_true',default=False,
        help='when used with -gau_temp, will return a gaussian input and exit')
    parser.add_argument('-func',     dest='functional',    type=str,  default='hf',help='specify the functional/method to be used, default=hf')
    parser.add_argument('-basis',    dest='basis',         type=str,  default='aug-cc-pvtz',help='specify the basis set to be used, default=aug-cc-pvtz')
    parser.add_argument('-noresp',   dest='noresp2',       action='store_true',default=True,help='skip the stage2 in resp py package, default=True')
    parser.add_argument('-scr',      dest='scratch',       type=str,  default='/dev/shm', help='specify the scratc dir for psi4, default=/dev/shm')
    parser.add_argument('-opdb',     dest='opdb',        type=str,  default='prefix',help='write the optimized coordinates, need to use together with -ipdb')
    parser.add_argument('-o',        dest='ofi',           type=str,  default='prefix',help='the psi4 output file, default = <input_prefix>.out')
    parser.add_argument('-odat',     dest='ome',           type=str,  default='prefix',help='file to record the results, default = <input_prefix>.dat')
    return parser

def read_xyz(xyz,chgmul='0,1',flag='psi4'):
    string = "{} \n".format(chgmul.replace(',', ' ').replace('e', '-'))
    pxyzs, gxyzs = [], ''
    units = False
    for line in open(xyz,'r'):
        #risky, if the comment line has more than 3 fields?
        if len(line.split()) > 3:
            string += line
            gxyzs += line
            pxyzs.append([float(x.strip()) for x in line.split()[1:]])
        if "units angstrom" in line:
            units = True
    if not units:
        string += "units angstrom\n"
    if flag == 'psi4': return string
    if flag == 'gau': return gxyzs
    if flag == 'pureXYZ': return pxyzs

def read_pdb(pdb):
    out, con = [], []
    for line in open(pdb,'r'):
        if line[0:6] in ['ATOM  ', 'HETATM']:
            out.append([line[0:6].strip(),int(line[6:11].strip()),line[12:16].strip().rstrip(),line[16:17],
                line[17:20].strip(),line[21:22],int(line[22:26].strip()),line[26:27],
                float(line[30:38].strip()),float(line[38:46].strip()),float(line[46:54].strip()),float(line[54:60].strip()),
                float(line[60:66].strip()),line[76:78].strip(),line[78:80]])
    #items in out: [0]atom/hetatm; [1]atom serial number; [2]atom name; 
    # [3]alternate location indicator; [4]residue name; [5]chain identifier
    # [6]residue sequence number; [7]code for insertion of residues (can be empty)
    # [8],[9],[10]orthogonal XYZ coords; [11]occupancy; [12] temperature factor;
    # [13]element symbol; [14]charge on the atom
        if 'CONECT' in line or 'END' in line: con.append(line)
    return [out, con]

def write_pdb(pdbinfo,coor,oname):
    #limited functions, just replacing the XYZs
    PDBForm = "{:6s}{:5d} {:^4s}{:1s}{:3s} {:1s}{:4d}{:1s}   {:8.3f}" + \
                    "{:8.3f}{:8.3f}{:6.2f}{:6.2f}          {:>2s}{:2s}\n"
    if len(pdbinfo) != len(coor):
        print("number of atom are not consistent!")
        print(len(pdbinfo), len(coor))
        exit(1)
    with open(oname,'w') as fp:
        for i in range(len(coor)):
            out = pdbinfo[i][0:8]+[float(x) for x in coor[i]]+pdbinfo[i][11:]
            fp.write(PDBForm.format(*out))
    return 0

def get_info(spe, wfn):
    Ehomo = wfn.epsilon_a_subset("AO","ALL").np[wfn.nalpha()-1]
    Elumo = wfn.epsilon_a_subset("AO","ALL").np[wfn.nalpha()]
    elecI = ((Ehomo+Elumo)/2)**2/(2*(Elumo-Ehomo))*27.211
    infor = f'''SCF(a.u),E_HOMO(a.u.),E_LUMO(a.u.),EI (eV)
{spe},{Ehomo},{Elumo},{elecI}
'''
    print(infor)
    return infor


# Main starts here:
args = read_par()
args = args.parse_args()
if args.ixyz != 'None' and args.ipdb == 'None':
    basename = '.'.join(args.ixyz.split('.')[0:-1])
    xyz = read_xyz(args.ixyz, args.chgmul)
    cor = read_xyz(args.ixyz,flag='gau')
    fname = '.'.join(args.ome.split('.')[0:-1])
elif args.ipdb != 'None' and args.ixyz == 'None':
    basename = '.'.join(args.ipdb.split('.')[0:-1])
    pdb = read_pdb(args.ipdb)[0]
    cor = '\n'.join(["{} {} {} {}\n".format(x[13],*x[8:11]) for x in pdb])
    xyz = "{}\n".format(args.chgmul.replace(',', ' ').replace('e', '-'))
    xyz += cor
    xyz += "units angstrom\n"
    if args.opdb == 'prefix': args.opdb = f'{basename}.opt.pdb'
elif args.ixyz != 'None' and args.ipdb != 'None':
    #just do the convert and exit
    xyz = read_xyz(args.ixyz,flag='pureXYZ')
    pdb = read_pdb(args.ipdb)[0]
    basename = '.'.join(args.ipdb.split('.')[0:-1])
    if args.opdb == 'prefix': args.opdb = basename+'.new.pdb'
    write_pdb(pdb,xyz,args.opdb)
    print('both -ipdb and -ixyz are given, I just write a new pdb with\n\
    name from -ipdb, coordinates from -ixyz')
    exit(0)
else:
    print("please provide -ixyz or/and -ipdb!")
    exit(1)

if args.gau:
    #writeh a template for gaussian opt and esp calculation
    gau_f = open("{}.com".format(fname), 'w')
    gau_i = (
    f"%chk={fname}.chk\n"
    f"%mem=64GB\n"
    f"%nproc=32\n"
    f"#p opt b3lyp/6-31g**\n\n"
    f"no titile\n\n"
    f"{args.chgmul.replace(',', ' ')}\n"
    f"{cor}\n"
    f"--link1--\n"
    f"%chk={fname}.chk\n"
    f"%nproc={str(args.core)}\n"
    f"%mem={str(args.memory)}GB\n"
    f"#p {args.functional}/{args.basis} pop=mk iop(6/33=2,6/50=1)\n"
    f"geom=check guess=read\n\n"
    f"no title\n\n"
    f"{args.chg} {args.mul}\n\n"
    f"{fname}.esp\n\n\n\n"
    )
    gau_f.write(gau_i)
    gau_f.close()

if args.dry:
    exit(0)

if args.ofi == 'prefix': args.ofi = f'{basename}.out'
if args.ome == 'prefix': args.ome = f'{basename}.dat'
psi4.core.set_output_file(args.ofi, True)
psi4_io = psi4.core.IOManager.shared_object()
psi4_io.set_default_path(args.scratch)

mol = psi4.geometry(xyz)
mol.update_geometry()
mol.set_name(basename)
out_f = open(args.ome, 'w')

psi4.set_memory(str(args.memory)+'GB')
psi4.set_num_threads(args.core)


if not args.no_opt:
    # quick optimization
    psi4.optimize("b3lyp/6-31g*", molecule=mol)
    mol.save_xyz_file(basename+"_opt.xyz",1)
    if args.opdb != 'None' and args.ipdb != 'None':
        coor = [x.split()[1:] for x in mol.save_string_xyz().split('\n') if len(x)>5]
        write_pdb(pdb,coor,args.opdb)
    scf_e, scf_wfn = psi4.energy("b3lyp/6-31g*", return_wfn=True)
    info = get_info(scf_e,scf_wfn)
    out_f.write('\nResults in gas phase:\n')
    out_f.write(info)
    mol.update_geometry()

## calculation of ESP:
options = {'VDW_SCALE_FACTORS'  : [1.4, 1.6, 1.8, 2.0],
           'VDW_POINT_DENSITY'  : 1.0,
           'RESP_A'             : 0.0005,
           'RESP_B'             : 0.1,
           'METHOD_ESP'         : args.functional,
           'BASIS_ESP'          : args.basis,
           }

# Call for first stage fit
charges1 = resp.resp([mol], options)
out_f.write('Electrostatic Potential Charges')
out_f.write('\n'.join([str(x) for x in charges1[0]]))
out_f.write('Restrained Electrostatic Potential Charges')
out_f.write('\n'.join([str(x) for x in charges1[1]]))
if args.noresp2:
    out_f.close()
    exit(0)

# Change the value of the RESP parameter A
# The second stage is problematic; because it change the ESP points of the first
# stage. So, default is to skip this step
options['RESP_A'] = 0.001

# Add constraint for atoms fixed in second stage fit
# some numbers need to be adjusted from antechamber
# And it may be useful in ligand_residue_cappers
#constraint_charge = []
#for i in range(4, 8):
#    constraint_charge.append([charges1[1][i], [i+1]])
#options['constraint_charge'] = constraint_charge
#options['constraint_group'] = [[2, 3, 4]]
options['grid'] = ['1_{:s}_grid.dat'.format(fname)]
options['esp'] = ['1_{:s}_grid_esp.dat'.format(fname)]

# Call for second stage fit
charges2 = resp.resp([mol], options)

# Get RESP charges
out_f.write("\nStage Two:\n")
out_f.write('RESP Charges')
out_f.write(charges2[1])
