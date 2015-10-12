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

def usage():
	print >> stderr, 'Usage: ', argv[0], '[--truncate-mprecision 0.95 | --ignore-class name_of_reject_class ] < matrices.cmat'


if __name__=="__main__":
	import getopt

	# parse command line options
	try:
		opts, args = getopt.getopt( argv[1:], 'hb:t:i:', ['help','truncate-mprecision=','ignore-class='] )
	except getopt.GetoptError, err:
		print str( err ) # will print something like "option -a not recognized"
		usage()
		exit(2)

	ignore_class = ""
	truncate = 0
		
	for o, a in opts:
		if o in ("-h", "--help"):
			usage()
			exit()
		elif o in ("-t", "--truncate-mprecision"):
			if "." in a:
				truncate = float( a )
			else:
				truncate = int( a )
		elif o in ("-i", "--ignore-class"):
			ignore_class = a
		else:
			assert False, "unhandled option"
	
	print "rank\tprecision\t precision std.\tnumber precision bins\trecall\trecall std.\tnumber real bins\taccuracy\tmisclassification rate"
	
	for cmat in parseConfusionMatrix( stdin ):
		
		acc = cmat.accuracy( ignore_class )*100
		mis = cmat.misclassification_rate( ignore_class )*100
		rec, rec_std, rec_num = cmat.macro_recall( ignore_class=ignore_class )
		prec, prec_std, prec_num = cmat.macro_precision( ignore_class=ignore_class, truncate=truncate )
		rec *= 100.
		rec_std *= 100.
		prec *= 100.
		prec_std *= 100.
		
		print "%s\t%.2f\t%.2f\t%i\t%.2f\t%.2f\t%i\t%.2f\t%.2f" % (cmat.title, prec, prec_std, prec_num, rec, rec_std, rec_num, acc, mis)
