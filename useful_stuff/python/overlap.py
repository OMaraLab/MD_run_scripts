#!/usr/bin/env python

import MDAnalysis as mda
import itertools
import argparse

"""
This is a python script to compare coordinate files and find resides present in
all structures. The intended use case was to make it easier to use vmd rmsd trajectory
tool to align multiple structures of the same protein, when these structures have small
sequence differences such as point mutations, or missing residues.  One example of this
would be aligning a cryoEM structure that is missing some residues in a flexible linker
with a full length alphafold structure prediction.

USE WITH CAUTION!  overlap.py is a very stupid program, and operates solely on residue
numbers.  It does not do any kind of sequence alignment, or consider residue similarity
in any way.

overlap.py does not consider ACE or NH2 caps to be part of a continuous protein sequence.

If you pass overlap.py a single coordinate file, it return a list of all contiguous 
sequences of amino acids in that file, and all discontinuities in the amino acid 
sequence, starting from the first protein residue in the file, and ending with the final
protein residue in the file. It will not identify missing residues at the termini.

Behaviour with coordinate files containing multiple sequences is not untested.

If you pass overlap.py multiple coordinate files, it will return a list of residue 
numbers which are present in the protein sequence of every coordinate file. 

overlap.py does not currently see point mutations as a mismatch. If residue 123 is MET 
in one structure but TYR in another, overlap.py considers the residue present in both.
"""

def loader(input_files):
    u = []
    p = []

    for i in input_files:
        if verbose_output:
            print(f'     Loading {i}')
        system = mda.Universe(i)
        u.append(system)
        p.append(system.select_atoms("protein and not resname ACE and not resname NH2"))
    if verbose_output:
        print('Done')

    return(p)

def to_ranges(iterable):
    iterable = sorted(set(iterable))
    # https://stackoverflow.com/questions/4628333/converting-a-list-of-integers-into-range-in-python
    for key, group in itertools.groupby(enumerate(iterable),
                                        lambda t: t[1] - t[0]):
        group = list(group)
        yield group[0][1], group[-1][1]

def vmdstr(a):
    selstr = "    "
    c = 0
    for i in list(to_ranges(a)):
        if c != 0: selstr += (" or ")
        selstr += (" ".join(["resid",str(i[0]),"to",str(i[1])]))
        c += 1
    return(selstr)

def overlap(p):
    reslist = [x.resnum for x in list(p[0].residues)]
    for i in p[1:]:
        a = [x.resnum for x in list(i.residues)]
        overlap = [x for x in a if x in reslist]
        reslist = overlap
    return (overlap)

def contiguous(p):
    a = [x.resnum for x in list(p[0].residues)]
    a
    r = range(min(a),max(a)+1)
    missing = [x for x in r if x not in a]
    present = [x for x in r if x in a]
    mlist = list(to_ranges(missing))
    plist = list(to_ranges(present))    
    return (missing, present)


# Create the argument parser object
parser = argparse.ArgumentParser(description="Find vmd strings to match sections of proteins.  If a single structure file is passed, find contiguous and missing regions in protein sequence.  If multiple coodinate files are passed, find protein residues shared between all. Only intended for files containing single polypeptide chains")

# Add arguments to the parser
parser.add_argument('input_files', type=str, nargs='+', help='Input file paths.  Each should contain a single polypeptide chain.')
parser.add_argument('--output', type=str, required=False, help='Output file path')
parser.add_argument('--verbose', action='store_true', help='Verbose output')

# Parse the arguments from the command line
args = parser.parse_args()

# Access the argument values
input_files = args.input_files
output_file = args.output
verbose_output = args.verbose

# Use the argument values in the script
if verbose_output:
    print(f'Finding VMD selection strings for {len(input_files)} input files.')

p = loader(input_files)

if len(p)==1:
    # if only passed one file, find contiguous sequences and discontinuities
    if verbose_output:
        print(f'\nOnly one file, identifying contiguous regions in protein structure.\n')
    missing, present = contiguous(p)
    vmd_p = vmdstr(present)
    vmd_m = vmdstr(missing)
    print('VMD selection string for all continuous blocks of residues:\n\n',vmd_p,'\n')
    print('VMD selection string for all discontinuities:\n\n',vmd_m)
    if output_file:
        with open(output_file, 'w') as f:
            f.write('\n'.join([f'Input files:']+["    " + x for x in input_files] 
                + ['',f'VMD selection string for all continuous blocks of residues:\n',vmd_p,''] 
                + [f'VMD selection string for all discontinuities in input files\n',vmd_m]
                ))

else:
    # if passed multiple files, find resids present in every file

    if verbose_output:
        print(f'\nExamining multiple structures to identify residues present in all structures.\n')
    olist = overlap(p)
    vmd_o = vmdstr(olist)
    print('VMD selection string for residues shared between all structures is:\n\n',vmd_o)
    if output_file:
        with open(output_file, 'w') as f:
            f.write('\n'.join([f'Input files:\n'] + ["    " + x for x in input_files]  + ["",'VMD selection string for overlapping residues:\n',vmd_o]))

if verbose_output:
    if output_file:
        print(f'\nSaving to {output_file}')
