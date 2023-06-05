#!/usr/bin/env python

import os
from argparse import ArgumentParser
from collections import defaultdict
from shutil import which
from sys import argv, exit, stderr, stdout

# Function definitions

def check_flag(flag):
    """Check the SAM flag and return booleans indicating the type of mapping"""
    flag = int(flag)
    proper_pair = (flag >> 1) & 1
    single = (flag >> 3) & 1
    not_primary = (flag >> 8) & 1
    return proper_pair, single, not_primary

# Define all command-line arguments

parser = ArgumentParser(description='Reads a sorted BAM file with "samtools -view -h bamfile" command and\
                                     a FASTA file of the reference sequences.  Options can set handling reads\
                                     mapping to multiple references, single vs paired reads, etc.  The\
                                     samtools command must be in the system path.  Output is tab-delimited,\
                                     two columns, reference name, TPM value.')

parser.add_argument('--fasta', metavar='', type=str,
                   help='FASTA file of reference sequences (required)')

parser.add_argument('--bam', metavar='', type=str,
                   help='Path to sorted BAM file (required)')

parser.add_argument('--value', metavar='', type=str,
                   help='Value to report (required; \'counts\' or \'tpm\')')

parser.add_argument('--output', metavar='', type=str, default=None,
                   help='Output file (optional, default=stdout)')

parser.add_argument('--decimals', metavar='', type=int, default=4,
                   help='Number of decimal places for TPM values (optional, default=4)')

parser.add_argument('--single', action='store_true', default=False,
                   help='Include single reads in TPM counts (optional, defaut=False)')

parser.add_argument('--progress', action='store_true', default=False,
                   help='Print progress to the screen (optional, default=False)')

# Print the help message if no arguments are provided, otherwise parse the arguments

if len(argv)==1:
    parser.print_help(stderr)
    exit(1)

args = parser.parse_args()

if args.value.lower() not in ('counts','tpm'):
    print('Error, --value must be either \'counts\' or \'tpm\' not %s; exiting.' % args.value.lower())
    exit(1)

# Global containers

lengths = defaultdict(int) # Dict of reference sequence lengths

for line in open(args.fasta):
    line = line.strip()
    if line and line[0] == '>':
        header = line[1:]
    else:
        lengths[header] += len(line)

og_counts = { x.split()[0]:0 for x in lengths.keys() } # Dict for read counts per ref sequence

done = set() # Set of concordantly mapped read pairs, to avoid double-counting these

# Check if samtools runs

if which('samtools') is None:
    print('Error: SamTools is not in the system path.  Exiting.')
    exit(1)

# Pipe the SamTools view command output and iterate over lines

if args.progress:
    print('Reading input...', file=stderr)

count = 0
    
with os.popen('samtools view -h %s' % args.bam) as pipe:
    for line in pipe:
        count += 1
        if args.progress:
            if count % 1000 == 0:
                print(' %sK lines processed\r' % (count // 1000), file=stderr, end='')
        if line[0] != '@':
            line = line.split()
            read, flag, ref, mate = line[0], line[1], line[2], line[6] 
            proper_pair, single, not_primary = check_flag(flag)
            if not_primary:
                continue
            if proper_pair and (mate=='=') and (read not in done):
                og_counts[ref] += 1
                done.add(read)
            elif single and args.single:
                og_counts[ref] += 1

if args.progress:
    print('Done.', file=stderr)          

# Set output to provided path or stdout

if args.output:
    out = open(args.output,'w')
else:
    out = stdout
    
# Print the counts if --value is counts and exit without doing TPM

if args.value.lower() == 'counts':
    if args.progress:
        print('Writing counts output...', file=stderr)
    
    for og,count in og_counts.items():
        print(og,count,sep='\t',file=out)
        
    if args.progress:
        print('Done.', file=stderr)
    
    exit()

# Otherwise, calculate and print TPM values
    
if args.progress:
    print('Calculating TPMs...', file=stderr)

og_rpks = { og:(count/lengths[og]) for og,count in og_counts.items() }
scaling = sum(og_rpks.values()) / 1e6
og_tpms = { og:(rpk/scaling) for og,rpk in og_rpks.items() }

if args.progress:
    print('Done.', file=stderr)

# Print the TPM values

if args.progress:
    print('Writing TPM output...', file=stderr)

for og,tpm in sorted(list(og_tpms.items())):
    print(og,round(tpm,args.decimals),sep='\t',file=out)
    
if args.progress:
    print('Done.', file=stderr)

