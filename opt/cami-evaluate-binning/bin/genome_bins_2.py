#!/usr/bin/python

import sys
import gzip
from exceptions import RuntimeError


def get_genome_mapping(mapping_file):
	"""
	# the mapping file needs to follow the exact format
	#anonymous_contig_id    genome_id       tax_id  contig_id       number_reads    start_position  end_position    total_length
	RL|S1|C10817    Sample18_57     45202   Sample18_57     44394   2       20519   20518

	@param mapping_file:
	@type mapping_file: str | unicode

	@return:
	"""
	with gzip.open(mapping_file) as mapping:
		genome_id_to_total_length = {}
		genome_id_to_list_of_contigs = {}
		anonymous_contig_id_to_genome_id = {}
		anonymous_contig_id_to_lengths = {}
		for line in mapping:
			if len(line.strip()) == 0 or line.startswith('#'):
				continue
			data = line.split('\t')
			if len(data) > 4:
				anonymous_contig_id, genome_id, tax_id, contig_id, number_reads, start_position, end_position, total_length = data
				total_length = float(total_length)
			else:
				anonymous_contig_id, genome_id, tax_id, contig_id = data
				total_length = 150. * 2  # doubled because mapping file contains paired reads

			anonymous_contig_id_to_lengths[anonymous_contig_id] = total_length
			anonymous_contig_id_to_genome_id[anonymous_contig_id] = genome_id
			if genome_id not in genome_id_to_total_length:
				genome_id_to_total_length[genome_id] = 0
				genome_id_to_list_of_contigs[genome_id] = []
			genome_id_to_total_length[genome_id] += total_length
			genome_id_to_list_of_contigs[genome_id].append(anonymous_contig_id)
		mapping.close()
	return genome_id_to_total_length, genome_id_to_list_of_contigs, anonymous_contig_id_to_genome_id, anonymous_contig_id_to_lengths


def read_header(input_stream):
	"""
	@Version:0.9.1
	@SampleID:RH_S1

	@@SEQUENCEID    BINID   TAXID
	"""
	header = {}
	column_names = {}
	for line in input_stream:
		if len(line.strip()) == 0 or line.startswith("#"):
			continue
		line = line.rstrip('\n')
		if line.startswith("@@"):
			for index, column_name in enumerate(line[2:].split('\t')):
				column_names[column_name] = index
			return header, column_names
		if line.startswith("@"):
			key, value = line[1:].split(':', 1)
			header[key] = value.strip()


def read_rows(input_stream, index_key, index_value):
	for line in input_stream:
		if len(line.strip()) == 0 or line.startswith("#"):
			continue
		line = line.rstrip('\n')
		row_data = line.split('\t')
		key = row_data[index_key]
		value = row_data[index_value]
		yield key, value


def read_participant_file(input_stream):
	"""
	# @@SEQUENCEID    BINID   TAXID

	@param: file_path
	@rtype: generator[dict[str|unicode, str|unicode]]
	"""
	header, column_names = read_header(input_stream)
	if "SEQUENCEID" not in column_names:
		raise RuntimeError("Column not found: {}".format("BINID"))

	if "BINID" in column_names:
		index_value = column_names["BINID"]
	elif "TAXID" in column_names:
		index_value = column_names["TAXID"]
	else:
		raise RuntimeError("Column not found: {}".format("BINID/TAXID"))

	index_key = column_names["SEQUENCEID"]

	return read_rows(input_stream, index_key, index_value)

def validate_genomes(file_path_query, file_path_mapping, file_path_output):
    """
    This script calculates precision, recall, accuracy and ARI for binned contigs
    
    @param file_path_query:
    @param file_path_mapping:
    @param file_path_output:
    @return:
    """
    #import ipdb; ipdb.set_trace()
    genome_id_to_total_length, genome_id_to_list_of_contigs, sequence_id_to_genome_id, anonymouse_contig_id_to_length = get_genome_mapping(file_path_mapping)
    sequence_id_to_bin_id = {}
    bin_id_to_list_of_sequence_id = {}
    bin_id_to_total_length = {}
    with open(file_path_query) as read_handler:
        for sequence_id, predicted_bin in read_participant_file(read_handler):
            sequence_id_to_bin_id[sequence_id] = predicted_bin
            if predicted_bin not in bin_id_to_total_length:
                bin_id_to_list_of_sequence_id[predicted_bin] = []
                bin_id_to_total_length[predicted_bin] = 0
            bin_id_to_list_of_sequence_id[predicted_bin].append(sequence_id)
            bin_id_to_total_length[predicted_bin] += anonymouse_contig_id_to_length[sequence_id]
    
    #now compute confusion matrix as dictionary of Counters    
    bin_id_to_genome_id_to_total_length = defaultdict(Counter)
    for predicted_bin in bin_id_to_list_of_sequence_id:
        for sequence_id in bin_id_to_list_of_sequence_id[predicted_bin]:
            genome_id = sequence_id_to_genome_id[sequence_id]
            bin_id_to_genome_id_to_total_length[predicted_bin][genome_id] += anonymouse_contig_id_to_length[sequence_id]

    (precision,recall) = calcPrecisionRecall(bin_id_to_genome_id_to_total_length, bin_id_to_total_length, genome_id_to_total_length)

    with open(file_path_output + "/by_genome_weighted.tsv", 'w') as write_handler:
        write_handler.write("Precision\tRecall\n")
        write_handler.write("%s \t %s\n" % (precision,recall))

def map_genomes(file_path_query, file_path_mapping, file_path_output):
	"""
	This script mapps a predicted bin to the genome with the highest recall

	@attention: In case of reads, read ids might not be paired read id and cause error: ReadID/1 ReadID/2

	@param file_path_query:
	@param file_path_mapping:
	@param file_path_output:
	@return:
	"""
	genome_id_to_total_length, genome_id_to_list_of_contigs, sequence_id_to_genome_id, anonymous_contig_id_to_lengths = get_genome_mapping(file_path_mapping)
	sequence_id_to_bin_id = {}
	bin_id_to_list_of_sequence_id = {}
	bin_id_to_total_lengths = {}
	with open(file_path_query) as read_handler:
		for sequence_id, predicted_bin in read_participant_file(read_handler):
			sequence_id_to_bin_id[sequence_id] = predicted_bin
			if predicted_bin not in bin_id_to_total_lengths:
				bin_id_to_list_of_sequence_id[predicted_bin] = []
				bin_id_to_total_lengths[predicted_bin] = 0
			bin_id_to_list_of_sequence_id[predicted_bin].append(sequence_id)
			bin_id_to_total_lengths[predicted_bin] += anonymous_contig_id_to_lengths[sequence_id]
	bin_metrics = {}
	bin_id_to_genome_id_to_total_length = {}
	mapped = set()
	for predicted_bin in bin_id_to_list_of_sequence_id:
		bin_id_to_genome_id_to_total_length[predicted_bin] = {}
		for sequence_id in bin_id_to_list_of_sequence_id[predicted_bin]:
			genome_id = sequence_id_to_genome_id[sequence_id]
			if genome_id not in bin_id_to_genome_id_to_total_length[predicted_bin]:
				bin_id_to_genome_id_to_total_length[predicted_bin][genome_id] = 0
			bin_id_to_genome_id_to_total_length[predicted_bin][genome_id] += anonymous_contig_id_to_lengths[sequence_id]
		max_overlap = 0.0
		best_genome_id = ""
		for genome_id in bin_id_to_genome_id_to_total_length[predicted_bin]:
			overlap = bin_id_to_genome_id_to_total_length[predicted_bin][genome_id]
			# bin_id_to_genome_id_to_total_length[predicted_bin][genome_id] /= genome_id_to_total_length[genome_id]
			if overlap > max_overlap:
				max_overlap = overlap
				best_genome_id = genome_id
		if best_genome_id not in mapped:
			mapped.add(best_genome_id)
		# length of genome in bin divided by bin size
		precision = bin_id_to_genome_id_to_total_length[predicted_bin][best_genome_id] / bin_id_to_total_lengths[predicted_bin]
		recall = bin_id_to_genome_id_to_total_length[predicted_bin][best_genome_id] / genome_id_to_total_length[best_genome_id]
		bin_metrics[predicted_bin] = [best_genome_id, precision, recall, genome_id_to_total_length[best_genome_id]]
	with open(file_path_output + "/by_genome.tsv", 'w') as write_handler:
		write_handler.write("Instance\tclass\tprecision\trecall\tpredicted class size\treal class size\n")
		for predicted_bin in bin_metrics:
			write_handler.write("strain\t%s\t%s\t%s\t%s\t%s\n" % (
				bin_metrics[predicted_bin][0],
				bin_metrics[predicted_bin][1],
				bin_metrics[predicted_bin][2],
				bin_id_to_total_lengths[predicted_bin],
				bin_metrics[predicted_bin][3]))
		for genid in genome_id_to_list_of_contigs:
			if genid not in mapped:
				write_handler.write("strain\t%s\t%s\t%s\t%s\t%s\n" % (genid, 'NA', 0., 0, genome_id_to_total_length[genid])) #precision is NA for unpredicted bins


map_genomes(
	file_path_query=sys.argv[1],
	file_path_mapping=sys.argv[2],
	file_path_output=sys.argv[3]
	)
	
validate_genomes(
	file_path_query=sys.argv[1],
	file_path_mapping=sys.argv[2],
	file_path_output=sys.argv[3]
	)
