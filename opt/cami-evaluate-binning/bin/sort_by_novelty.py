#!/usr/bin/env python

import sys
import os
import tempfile
import argparse


def read_novelties(path_novelty):
	gid_to_novelty = {}
	with open(path_novelty) as read_handler:
		for line in read_handler:
			line = line.rstrip('\n')
			if line.startswith('#') or len(line) == 0:
				continue
			gid, novelty = line.split('\t')
			gid_to_novelty[gid] = novelty
	return gid_to_novelty


def parse_cami_format(file_path, required_columns, key_column=None):
	"""

	@type file_path: str
	@type required_columns: list[str]
	@type key_column: str
	@rtype: dict[str, dict[str, str]]
	"""
	if key_column is None:
		key_column = "SEQUENCEID"
	with open(file_path) as read_handler:
		column_names = []
		for line in read_handler:
			line = line.strip()
			if line.startswith('@@'):
				column_names = line[2:].split('\t')
				break
			if line.startswith('@') or line.startswith('#') or len(line) == 0:
				continue

		assert len(column_names) > 0, "No column declaration found."
		for column in required_columns:
			assert column in column_names, "Column '{}' not found in '{}'.".format(column, os.path.basename(file_path))
		key_index = column_names.index(key_column)

		for line in read_handler:
			line = line.rstrip('\n')
			if line.startswith('#') or len(line) == 0:
				continue
			row = line.split('\t')
			# cami_row_dict[row[key_index]] = {}
			cami_row_dict = {}
			for column in required_columns:
				value_index = column_names.index(column)
				cami_row_dict[column] = row[value_index]
			yield row[key_index], cami_row_dict
# 	return cami_dict


def split_seq_file(sid_to_novelty, file_path_input, out_dir):
	dict_file_handlers = {}
	dict_file_paths = {}
	with open(file_path_input) as read_handler:
		for line_seq_length in read_handler:
			sid, _ = line_seq_length.split('\t')
			novelty = sid_to_novelty[sid]
			if novelty not in dict_file_handlers:
				file_path = tempfile.mktemp(prefix="seq_len_", dir=out_dir)
				dict_file_paths[novelty] = file_path
				dict_file_handlers[novelty] = open(file_path, 'w')
			dict_file_handlers[novelty].write(line_seq_length)

	for novelty in dict_file_handlers:
		dict_file_handlers[novelty].close()
	return dict_file_paths


def split_gold_file(sid_to_novelty, file_path_input, out_dir):
	cami_format_header_gold = "@Version:0.9.1\n\n@@SEQUENCEID\tTAXID\tBINID\n"
	dict_file_handlers = {}
	dict_file_paths = {}
	for sid, columns in parse_cami_format(file_path_input, ["SEQUENCEID", "TAXID", "BINID"], "SEQUENCEID"):
		bid = columns["BINID"]
		tid = columns["TAXID"]
		novelty = sid_to_novelty[sid]
		if novelty not in dict_file_handlers:
			file_path_gold = tempfile.mktemp(prefix="gold_", dir=out_dir)
			dict_file_paths[novelty] = file_path_gold
			dict_file_handlers[novelty] = open(file_path_gold, 'w')
			dict_file_handlers[novelty].write(cami_format_header_gold)
		dict_file_handlers[novelty].write("{sid}\t{tid}\t{bid}\n".format(sid=sid, tid=tid, bid=bid))

	for novelty in dict_file_handlers:
		dict_file_handlers[novelty].close()
	return dict_file_paths


def split_query_file(sid_to_novelty, file_path_input, out_dir, is_supervised, novelties):
	bin_column = "BINID"
	if is_supervised:
		bin_column = "TAXID"
	cami_format_header_query = "@Version:0.9.1\n\n@@SEQUENCEID\t{}\n".format(bin_column)
	dict_file_handlers = {}
	dict_file_paths = {}

	# making sure all novelty files will exist:
	for novelty in novelties:
		file_path_query = tempfile.mktemp(prefix="query_", dir=out_dir)
		dict_file_paths[novelty] = file_path_query
		dict_file_handlers[novelty] = open(file_path_query, 'w')
		dict_file_handlers[novelty].write(cami_format_header_query)

	for sid, columns in parse_cami_format(file_path_input, ["SEQUENCEID", bin_column], "SEQUENCEID"):
		novelty = sid_to_novelty[sid]
		dict_file_handlers[novelty].write("{sid}\t{bid}\n".format(sid=sid, bid=columns[bin_column]))

	for novelty in dict_file_handlers:
		dict_file_handlers[novelty].close()
	return dict_file_paths


def main(path_novelty, path_gold, path_query, path_seq_length, is_supervised, out_dir):
	assert path_novelty is not None, "Need novelty file"
	assert path_gold is not None, "Need gold standard file"
	assert path_query is not None, "Need query file"
	assert path_seq_length is not None, "Need seq length file"
	assert out_dir is not None, "Need output dir"

	sid_to_novelty = {}
	gid_to_novelty = read_novelties(path_novelty)
	for sid, columns in parse_cami_format(path_gold, ["SEQUENCEID", "TAXID", "BINID"], "SEQUENCEID"):
		sid_to_novelty[sid] = gid_to_novelty[columns["BINID"]]
	del gid_to_novelty

	novelty_to_file_paths_seq_length = split_seq_file(sid_to_novelty, path_seq_length, out_dir)
	novelty_to_file_paths_gold = split_gold_file(sid_to_novelty, path_gold, out_dir)
	novelties = novelty_to_file_paths_gold.keys()
	novelty_to_file_paths_query = split_query_file(sid_to_novelty, path_query, out_dir, is_supervised, novelties)

	for novelty in novelties:
		sys.stdout.write("{novelty}\t{gold}\t{query}\t{seql}\n".format(
			novelty=novelty,
			gold=novelty_to_file_paths_gold[novelty],
			query=novelty_to_file_paths_query[novelty],
			seql=novelty_to_file_paths_seq_length[novelty])
		)


def get_parser_options(args=None, version="Prototype"):
	"""
	Parsing of passed arguments.

	@param args: Passed arguments

	@return: any
	"""
	parser = argparse.ArgumentParser(
		usage="python %(prog)s configuration_file_path",
		version="MetagenomeSimulationPipeline {}".format(version),
		description="""""".format(),
		formatter_class=argparse.RawTextHelpFormatter)

	parser.add_argument(
		"-tax", "--supervised",
		action='store_true',
		default=False,
		help="")

	parser.add_argument(
		"-n", "--path_novelty",
		default=None,
		type=str,
		help="file path")

	parser.add_argument(
		"-g", "--path_gold",
		default=None,
		type=str,
		help="file path")

	parser.add_argument(
		"-q", "--path_query",
		default=None,
		type=str,
		help="file path")

	parser.add_argument(
		"-s", "--path_seq_length",
		default=None,
		type=str,
		help="file path")

	parser.add_argument(
		"-o", "--path_dir_out",
		default=None,
		type=str,
		help="file path")

	if args is None:
		return parser.parse_args()
	else:
		return parser.parse_args(args)


if __name__ == "__main__":
	options = get_parser_options()
	main(
		options.path_novelty, options.path_gold, options.path_query,
		options.path_seq_length, options.supervised, options.path_dir_out)
