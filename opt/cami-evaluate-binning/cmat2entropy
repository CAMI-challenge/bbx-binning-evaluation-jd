#!/usr/bin/env python
# -*- coding: utf-8 -*-
# This script takes confusion matrices from the standard input
# in the form:
#
# title\tcol1name\tcol2name\tcol3name
# row1name\tfreq11\tfreq12\freq13
# row2name\tfreq21\tfreq22\freq23
#
# where \t means TAB separation and there can be as many rows/columns as you
# like.
#
# From the matrices it will compute common statistics to measure classification
# accuracy and output them in the input file order.
#
# Conventions:
# a) First upper left cell can hold a title (any string)
# b) There is a special class name that will be ignored in some of the
#    statistics (defaults to empty string "")
# c) Any empty line marks the end of the last matrix and following lines are
#    interpreted as a new matrix.
# d) Lines starting with "#" are considered comments and are allowed anywhere
#    but inside the matrices 

from sys import argv, stdout, stderr, stdin, exit
from classevaltools import parseConfusionMatrix
from math import ceil


def usage():
    print >> stderr, 'Usage: ', argv[0], '[--ignore-class name_of_reject_class ] < matrices.cmat'


if __name__ == "__main__":
    import getopt

    # parse command line options
    try:
        opts, args = getopt.getopt(argv[1:], 'hi:', ['help', 'ignore-class='])
    except getopt.GetoptError, err:
        print str(err) # will print something like "option -a not recognized"
        usage()
        exit(2)

    ignore_class = ""
    truncate = 1.

    for o, a in opts:
        if o in ("-h", "--help"):
            usage()
            exit()
        elif o in ("-i", "--ignore-class"):
            ignore_class = a
        else:
            assert False, "unhandled option"

    print "instance\ttotal entropy"

    for cmat in parseConfusionMatrix(stdin):
        print "%s\t%.2f" % (cmat.title, cmat.entropy())
