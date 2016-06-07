#!/usr/bin/python

import sys

"""
Given a binning result (CAMI format), calculates a corresponding profile
Use in biobox, apart from that input is supposed to be in sys.argv:
[0]: name of this script
[1]: name of participant (filepath + name)
[2]: name of goldstandard (filepath + name)
[3]: Path to NCBI taxdump folder
[4]: output file path + name
[5]: contigbased: True/readbased: False
[6]: taxonomic: True/unsupervised: False
"""

# global variable
taxonomy = ['species','genus','family','order','class','phylum','superkingdom']

def map_genomes(filename, gsa):
	participant = open(filename)
	gold_standard = get_mapping(gsa, False)[2]
	start = False
	bins = {}
	seqs = {}
	for line in participant:
		if line.startswith('@@'):
			start = True
		elif start:
			line = line.split()
			if len(line) == 0:
				continue
			seqid = line[0] # Sequence ID
			if seqid in gold_standard:
				taxon = gold_standard[seqid]
			else: #shouldnt happen
				continue
			if line[1] in bins:
				seqs[line[1]].append(seqid)
				if taxon[0] in bins[line[1]]:
					bins[line[1]][taxon[0]] += taxon[1]
				else:	
					bins[line[1]][taxon[0]] = taxon[1]
			else:
				seqs[line[1]] = [line[0]]
				bins[line[1]] = {taxon[0] : taxon[1]}
	participant.close()
	gen_to_taxid = {}
	seq_to_taxid = {}
	for gbin in bins:
		max_len = 0
		max_tax = ""
		for val in bins[gbin]:
			if bins[gbin][val] > max_len:
				max_len = bins[gbin][val]
				max_val = val
		gen_to_taxid[gbin] = max_val
		seq_to_taxid[gbin] = max_val
	return (seqs,gen_to_taxid,gen_to_taxid)

def write_taxonomic(filename, gsa):
	bins,gen,seq_to_taxid = map_genomes(filename,gsa)
	p = open(filename)
	toWrite = ""
	for line in p: #write header
		toWrite += line
		if line.startswith('@@'):
			break
	for gbin in bins: #map seqid to taxid
		for seq in bins[gbin]:
			toWrite += seq
			toWrite += '\t'
			toWrite += gen[gbin]
			toWrite += '\n'
	return toWrite

profile = write_taxonomic(sys.argv[1],sys.argv[2])
name = sys.argv[1].split("/")[-1]
f = open("%s_t" % name,'w') # not beautiful but effective
f.write(profile)
	
