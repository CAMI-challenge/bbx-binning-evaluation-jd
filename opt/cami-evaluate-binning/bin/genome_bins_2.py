#!/usr/bin/python

import sys

def get_genomes(mapping_file): # the mapping file needs to follow the exact format
	mapping = open(mapping_file)
	gen_length = {}
	gen_map = {}
	c_gen = {}
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
		c_gen[contig] = gen_id
		if gen_id in gen_length:
			gen_length[gen_id] += length
			gen_map[gen_id].append(contig)
		else:
			gen_length[gen_id] = length
			gen_map[gen_id] = [contig]
	mapping.close()
	return (gen_length, gen_map, c_gen, lengths)

def map_genomes(filename, gsa, out):
	participant = open(filename)
	gen_length, gen_map, c_gen, lengths = get_genomes(gsa)
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
				seqs[binid].append(seqid)
				bin_lengths[binid] += lengths[seqid]
			else:
				seqs[binid] = [seqid]
				bin_lengths[binid] = lengths[seqid]
	participant.close()
	bin_metrics = {}
	bin_gen = {}
	mapped = []
	for binid in seqs:
		bin_gen[binid] = {}
		for seq in seqs[binid]:
			gen = c_gen[seq]
			if gen in bin_gen[binid]:
				bin_gen[binid][gen] += lengths[seq]
			else:
				bin_gen[binid][gen] = lengths[seq]
		max_rec = 0.0
		mgen = ""
		for gen in bin_gen[binid]:
			bin_gen[binid][gen] /= gen_length[gen]
			if bin_gen[binid][gen] > max_rec:
				max_rec = bin_gen[binid][gen]
				mgen = gen
		if mgen not in mapped:
			mapped.append(mgen)
		prec = bin_gen[binid][mgen] * gen_length[mgen] / bin_lengths[binid] # length of genome divided by bin size
		bin_metrics[binid] = [mgen,prec,max_rec,gen_length[mgen]]
	f = open(out,'w')
	toWrite = "Instance\tclass\tprecision\trecall\tpredicted class size\treal class size\n"
	for binid in bin_metrics:
		toWrite += "strain\t%s\t%s\t%s\t%s\t%s\n" % (bin_metrics[binid][0],bin_metrics[binid][1], bin_metrics[binid][2],bin_lengths[binid],bin_metrics[binid][3])
	for genid in gen_map:
		if genid not in mapped:
			toWrite += "strain\t%s\t%s\t%s\t%s\t%s\n" % (genid,0,0,0,gen_length[genid])
	f.write(toWrite)
	f.close()


map_genomes(sys.argv[1],sys.argv[2],sys.argv[3])
