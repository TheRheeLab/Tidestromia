#!/usr/bin/env python

from numpy import mean
from sys import argv, stderr, exit

if len(argv) == 1:
    print('\nThis script accepts a BLAST tabular output file containing all pairwise comparisons between two species,')
    print('excluding self-blast (species against itself).  BLAST output files should be concatenated together into a')
    print('single file.  Output is printed to standard output, and is tab delimited, with the following fields:\n')
    print('Header1, Header2, Avg.E-Value, Avg.Bit Score\n')
    print('USAGE: ReciprocalBestHits.py BlastFile\n')
    exit()

fname = argv[1]

queries_to_best_matches = {}

with open(fname,'r') as file:
    for line in file:
        line = line.split()
        query, match = line[:2]
        ev, bs = float(line[10]), float(line[11])
        if query not in queries_to_best_matches:
            queries_to_best_matches[query] = (match,ev,bs)

done = set()

for query, (match,ev,bs) in queries_to_best_matches.items():
    if match in done:
        continue
    if match in queries_to_best_matches:
        rcp_match, rcp_ev, rcp_bs = queries_to_best_matches[match]
    if rcp_match == query:
        ev, bs = mean((ev,rcp_ev)), mean((ev,rcp_bs))
        print(query,match,ev,bs,sep='\t')
        printed.update((query,match))

