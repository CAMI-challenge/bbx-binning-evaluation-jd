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


def read_cami_format(file_path, required_columns, key_column=None):
	"""

	@type file_path: str
	@type required_columns: list[str]
	@type key_column: str
	@rtype: dict[str, dict[str, str]]
	"""
	if key_column is None:
		key_column = "SEQUENCEID"
	cami_dict = {}
	with open(file_path) as read_handler:
		column_names = []
		for line in read_handler:
			line = line.rstrip('\n')
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
			cami_dict[row[key_index]] = {}
			for column in required_columns:
				value_index = column_names.index(column)
				cami_dict[row[key_index]][column] = row[value_index]
	return cami_dict


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
		"-o", "--path_dir_out",
		default=None,
		type=str,
		help="file path")

	if args is None:
		return parser.parse_args()
	else:
		return parser.parse_args(args)


def main(path_novelty, path_gold, path_query, is_supervised, out_dir):
	assert path_novelty is not None, "Need novelty file"
	assert path_gold is not None, "Need gold standard file"
	assert path_query is not None, "Need query file"
	assert out_dir is not None, "Need output dir"

	query_bin_column = "BINID"
	if is_supervised:
		query_bin_column = "TAXID"
	gid_to_novelty = read_novelties(path_novelty)
	cami_dict_gold = read_cami_format(path_gold, ["SEQUENCEID", "TAXID", "BINID"], "SEQUENCEID")
	cami_dict_query = read_cami_format(path_query, ["SEQUENCEID", query_bin_column], "SEQUENCEID")

	cami_format_header_query = "@Version:0.9.1\n\n@@SEQUENCEID\t{}\n".format(query_bin_column)
	cami_format_header_gold = "@Version:0.9.1\n\n@@SEQUENCEID\tTAXID\tBINID\n"
	novelty_to_out_file_handler_query = {}
	novelty_to_out_file_paths_query = {}
	novelty_to_out_file_handler_gold = {}
	novelty_to_out_file_paths_gold = {}
	for sid in cami_dict_gold:
		bid_gold = cami_dict_gold[sid]["BINID"]
		tid_gold = cami_dict_gold[sid]["TAXID"]
		bin_query = None
		if sid in cami_dict_query:
			bin_query = cami_dict_query[sid][query_bin_column]
		novelty = gid_to_novelty[bid_gold]
		if novelty not in novelty_to_out_file_handler_query:
			file_path_query = tempfile.mktemp(prefix="query_", dir=out_dir)
			file_path_gold = tempfile.mktemp(prefix="gold_", dir=out_dir)
			novelty_to_out_file_paths_query[novelty] = file_path_query
			novelty_to_out_file_paths_gold[novelty] = file_path_gold
			novelty_to_out_file_handler_query[novelty] = open(file_path_query, 'w')
			novelty_to_out_file_handler_gold[novelty] = open(file_path_gold, 'w')
			novelty_to_out_file_handler_query[novelty].write(cami_format_header_query)
			novelty_to_out_file_handler_gold[novelty].write(cami_format_header_gold)
		if bin_query is not None:
			novelty_to_out_file_handler_query[novelty].write("{sid}\t{tid}\n".format(sid=sid, tid=bin_query))
		novelty_to_out_file_handler_gold[novelty].write("{sid}\t{tid}\t{bid}\n".format(sid=sid, tid=tid_gold, bid=bid_gold))

	for novelty in novelty_to_out_file_handler_query:
		novelty_to_out_file_handler_query[novelty].close()
		novelty_to_out_file_handler_gold[novelty].close()
		sys.stdout.write("{novelty}\t{query}\t{gold}\n".format(
			novelty=novelty,
			query=novelty_to_out_file_paths_query[novelty],
			gold=novelty_to_out_file_paths_gold[novelty])
		)


if __name__ == "__main__":
	options = get_parser_options()
	main(options.path_novelty, options.path_gold, options.path_query, options.supervised, options.path_dir_out)

