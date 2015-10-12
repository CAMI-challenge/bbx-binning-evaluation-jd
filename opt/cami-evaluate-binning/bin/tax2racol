#!/usr/bin/env python
# This script generates a tab-separated file with taxonomic path information
# about predicted sequences or other entities. The first column contains the
# name of the sequence/entity and each column starting from column two holds
# the name of consecutively higher ranking taxa (NCBI taxonomy). Each column
# corresponds to a common rank, empty fields mean that there is no taxon
# assigned at the given rank.
#
# Example:
# entity_name\trank_name_A\trank_name_B...
#
# where \t means TAB separation and after the entity name there can be as many
# columns as you which but at least one.
#
# Input:
# The input is analogous to the output but restricted to only two columns,
# the first being the sequence identifier and the second an NCBI taxonomic ID.
# This corresponds to common output of taxonomic assignment tools like MEGAN or
# PhyloPythia(S)
#
# Conventions:
# a) Comment lines in input must start with '#' (first character) 
# b) The first output line, if starting with '#' will give the name of the ranks

# suppress warnings with TaxonomyNcbi package
import warnings
with warnings.catch_warnings():
    warnings.filterwarnings("ignore", category=DeprecationWarning)
    from taxonomyncbi import TaxonomyNcbi


def pathIterID(taxonomy, taxid):
    while taxid:
        rank = taxonomy.getRank(taxid)
        yield str(taxid), rank
        taxid = taxonomy.getParentNcbid(taxid)


def pathIterName(taxonomy, taxid):
    while taxid:
        name = taxonomy.getScientificName(taxid)
        rank = taxonomy.getRank(taxid)
        yield name, rank
        taxid = taxonomy.getParentNcbid(taxid)


def taxID2Ranks(taxonomy, rank2pos, pathIter, taxid): #TODO: alternative is dictionary
    path = [""] * len(rank2pos)
    for name, rank in pathIter(taxonomy, taxid):
        try:
            i = rank2pos[rank]
            if path[i]:
                stderr.write("Warning: multiple '%s' ranks in path: replace '%s' -> '%s'.\n" % (rank, path[i], name))
            path[i] = name
        except KeyError:
            pass
    return path


header = lambda ranks: "#identifier\t%s" % ("\t".join(ranks))


# print information on usage
def Usage():
    print >> stderr, 'Usage: ', argv[0], '--taxonomy-backend-file ncbitax_sqlite.db --ranks genus,family,order'


if __name__ == "__main__":
    from sys import stdin, stderr, exit, argv
    import getopt

    # parse command line options
    try:
        opts, args = getopt.getopt(argv[1:], 'ht:r:i', ['help', 'taxonomy-backend-file=', 'ranks=', 'show-identifiers'])
    except getopt.GetoptError, err:
        # print help information and exit:
        print str(err)  # will print something like "option -a not recognized"
        Usage()
        exit(2)

    # defaults
    ranks = ["species", "genus", "family", "order", "class", "phylum", "superkingdom"]
    path_iterator = pathIterName
    taxonomy_filename = None

    # option parsing
    for o, a in opts:
        if o in ("-h", "--help"):
            Usage()
            exit()
        elif o in ("-t", "--taxonomy-backend-file"):
            taxonomy_filename = a
        elif o in ("-r", "--ranks"):
            ranks = a.split(",")
        elif o in ("-i", "--show-identifiers"):
            path_iterator = pathIterID
        else:
            assert False, "unhandled option"

    if not taxonomy_filename:
        print >> stderr, "you must specify a taxonomy file"
        Usage()
        exit(3)

    # handle broken pipes
    import signal

    signal.signal(signal.SIGPIPE, signal.SIG_DFL)

    # be verbose
    print >> stderr, 'Using taxonomy file %s' % taxonomy_filename
    print >> stderr, 'Using ranks %s' % ",".join(ranks)

    tax = TaxonomyNcbi(taxonomy_filename, ranks)
    rank2pos = dict((v, i) for i, v in enumerate(ranks))

    print header(ranks)

    for line in stdin:
        if line[0] != "#":
            line = line.rstrip()
            try:
                ident, taxid = line.split("\t")[:2]
            except ValueError:
                stderr.write("error parsing, skipping line \"%s\"" % line)
            path = taxID2Ranks(tax, rank2pos, path_iterator, taxid)
            print "%s\t%s" % (ident, "\t".join(path))
