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
    print >> stderr, 'Usage: ', argv[
        0], '[--truncate-mprecision 0.95 | --ignore-class name_of_reject_class ] < matrices.cmat'


if __name__ == "__main__":
    import getopt

    # parse command line options
    try:
        opts, args = getopt.getopt(argv[1:], 'hb:t:i:', ['help', 'truncate-mprecision=', 'ignore-class='])
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
        elif o in ("-t", "--truncate-mprecision"):
            if "." in a:
                truncate = float(a)
            else:
                truncate = int(a)
        elif o in ("-i", "--ignore-class"):
            ignore_class = a
        else:
            assert False, "unhandled option"

    print "instance\tclass\tprecision\trecall\tpredicted class size\treal class size"

    upsize = ursize = 0
    for cmat in parseConfusionMatrix(stdin):
        totalsize = 0
        tmpstore = []
        for name, psize, pcorrect in cmat.precision_freqs():
            if name != ignore_class:
                prec = pcorrect / float(psize)
                rsize, rcorrect = cmat.recall_freq(name)
                assert rcorrect == pcorrect
                rec = rcorrect / float(rsize)
                totalsize += psize
                tmpstore.append((psize, rsize, prec, rec, name))
            else:
                upsize = psize
                ursize = cmat.recall_freq(name)[0]
        tmpstore.sort(reverse=True) #reverse sort from high to low bins

        if type(truncate) == float and truncate <= 1.:
            lastsize = cumsize = 0
            threshold = ceil(totalsize * truncate)
            noiseclass = False
            for psize, rsize, prec, rec, name in tmpstore:
                if noiseclass:
                    upsize += psize
                    ursize += rsize
                    continue
                if cumsize > threshold and psize < lastsize: ##treat equal size classes
                    upsize += psize
                    ursize += rsize
                    noiseclass = True
                    continue
                print "%s\t%s\t%.2f\t%.2f\t%i\t%i" % (cmat.title, name, prec, rec, psize, rsize)
                cumsize += psize
                lastsize = psize
            if upsize or ursize:
               print "%s\t%s\t%.2f\t%.2f\t%i\t%i" % (cmat.title, ignore_class, float("nan"), float("nan"), upsize, ursize)
        else:
            raise TypeError("truncate must either be a valid fraction (float) between 0 and 1")
