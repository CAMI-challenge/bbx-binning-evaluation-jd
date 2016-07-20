#!/usr/bin/python

import sys

def get_genomes(mapping_file): # the mapping file needs to follow the exact format
	mapping = open(mapping_file)
	gen_length = {}
	gen_map = {}
	lengths = {}
	for line in mapping:
		if (line.startswith('#')):
			continue
		parts = line.split()
		contig = parts[0]
		gen_id = parts[1]
		try:
			length = float(parts[3].split('_')[-1]) # name has _contiglength as last component
		except:
			length = 150. # read based
		lengths[contig] = length
		if gen_id in gen_length:
			gen_length[gen_id] += length
			gen_map[gen_id].append(contig)
		else:
			gen_length[gen_id] = length
			gen_map[gen_id] = [contig]
	mapping.close()
	return (gen_length, gen_map, lengths)

def map_genomes(filename, gsa, out):
	participant = open(filename)
	gen_length, gen_map, lengths = get_genomes(gsa)
	start = False
	bins = {}
	seqs = {}
	bin_lengths = {}
	for line in participant:
		if line.startswith('@@'):
			start = True
		elif start:
			line = line.split()
			if len(line) == 0:
				continue
			seqid = line[0] # Sequence ID
			binid = line[1] # predicted bin
			bins[seqid] = binid
			if binid in bin_lengths:
				bin_lengths[binid] += lengths[seqid]
			else:
				bin_lengths[binid] = lengths[seqid]
	participant.close()
	bin_metrics = {}
	for gen in gen_map:
		mapped = {}
		unmapped = 0
		unmapped_seqs = []
		for contig in gen_map[gen]:
			try:
				binid = bins[contig]
			except KeyError:
				unmapped += lengths[contig]
				unmapped_seqs.append(contig)
				continue
			if binid in mapped:
				mapped[binid] += lengths[contig]
			else:
				mapped[binid] = lengths[contig]
		if len(mapped) == 0:
			bin_metrics[gen] = (0,0,"Unassigned")
			continue
		max_recall = 0.
		mapped_bin = ""
		for binid in mapped:
			recall = mapped[binid]/gen_length[gen]
			if recall > max_recall:
				max_recall = recall
				mapped_bin = binid
		bin_metrics[gen] = (max_recall,mapped[mapped_bin]/bin_lengths[mapped_bin],mapped_bin)
	f = open(out,'w')
	toWrite = "Instance\tclass\tprecision\trecall\tpredicted class size\treal class size\n"
	for gen in bin_metrics:
		try:
			toWrite += "strain\t%s\t%s\t%s\t%s\t%s\n" % (gen,bin_metrics[gen][1], bin_metrics[gen][0],bin_lengths[bin_metrics[gen][2]],gen_length[gen])
		except KeyError:
			toWrite += "strain\t%s\t%s\t%s\t%s\t%s\n" % (gen,0,0,0,gen_length[gen])
	f.write(toWrite)
	f.close()


map_genomes(sys.argv[1],sys.argv[2],sys.argv[3])
