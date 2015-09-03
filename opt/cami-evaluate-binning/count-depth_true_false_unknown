#!/usr/bin/env python

from sys import argv, stdout, stderr, exit
from numpy import mean

# simple dummy weight function counting each sequences as one
class oneweight:
	__getitem__ = lambda self,key: 1

def usage():
	print >> stderr, 'Usage: ', argv[0], '--labels lab.racol --predictions pred.racol [--with-unknown-labels  --weights sequences.weights --scale .001]'

if __name__ == "__main__":
	import getopt
	
	# parse command line options
	try:
		opts, args = getopt.getopt( argv[1:], 'hl:p:w:s:u', ['help', 'labels=','predictions=','weights=','scale=','with-unknown-labels'] )
	except getopt.GetoptError, err:
		print str( err ) # will print something like "option -a not recognized"
		usage()
		exit(2)
	
	# default parameters
	reffile = None
	predfile = None
	weightfile = None
	unknown_labels = False
	scale = 1
	
	for o, a in opts:
		if o in ("-h", "--help"):
			usage()
			exit()
		elif o in ("-l", "--labels"):
			reffile = a
		elif o in ("-p", "--predictions"):
			predfile = a
		elif o in ("-w", "--weights"):
			weightfile = a
		elif o in ("-u", "--with-unknown-labels"):
			unknown_labels = True
		elif o in ("-s", "--scale"):
			scale = float( a )
		else:
			assert False, "unhandled option"
	
	if not reffile:
		print >>stderr, "you must specify a file for taxonomic labels"
		usage()
		exit( 3 )
	
	if not predfile:
		print >>stderr, "you must specify a file for taxonomic predictions"
		usage()
		exit( 4 )
	
	# read ref assignments
	ref={}
	with open( reffile, "r" ) as f:
		for line in f:
			if line[0] != "#":
				line = line.rstrip( "\n" ).split( "\t" )
				ref[line[0]] = line[1:]

	# read predictions
	pred={}
	with open( predfile, "r" ) as f:
		for line in f:
			if line[0] != "#":
				line = line.rstrip( "\n" ).split( "\t" )
				pred[line[0]] = line[1:]
	
	# read weights if given
	if weightfile:
		weight = {}
		with open( weightfile, "r" ) as f:
			for line in f:
				name, w = line.strip().split( "\t", 2 )[:2]
				weight[name] = int( w )
	else:
		weight = oneweight()
	
	# output only false lines in modified format
	correct = {}
	incorrect = {}
	unknowns = {}
	depth = 0
	counter = 0

	for seq, path in ref.items():
		try:
			counter += 1
			#print path, pred[seq]
			plen = min( len( path ), len( pred[seq]) )
			for otax, ptax in zip( path, pred[seq] ):
				if ptax == "":
					plen -= 1
				elif unknown_labels and otax == "":
					try:
						unknowns[plen] += weight[seq]
					except KeyError:
						unknowns[plen] = weight[seq]
					break
				elif ptax == otax:
					try:
						correct[plen] += weight[seq]
					except KeyError:
						correct[plen] = weight[seq]
					break
				else:
					try:
						incorrect[plen] += weight[seq]
					except KeyError:
						incorrect[plen] = weight[seq]
					break
			if not plen:
				try:
					correct[plen] += weight[seq]
				except KeyError:
					correct[plen] = weight[seq]
			depth = max( depth, plen )
			del pred[seq] #remove entry from predictions
		except KeyError: #handle no prediction as root assignment
			try:
				correct[0] += weight[seq]
			except KeyError:
				correct[0] = weight[seq]
			stderr.write( "%s not found in prediction file\n" % (seq) )

	for seq, path in pred.items():
		counter += 1
		plen = len( path )
		for tax in path:
			if tax != "":
				try:
					unknowns[plen] += weight[seq]
				except KeyError:
					unknowns[plen] = weight[seq]
				break
			else:
				plen -= 1
		if not plen:
			try:
				unknowns[0] += weight[seq]
			except KeyError:
				unknowns[0] = weight[seq]
		depth = max( depth, plen )

	if type( weight ) == dict:
		assert sum( correct.values() + incorrect.values() + unknowns.values() ) == sum( weight.values() )
	else:
		assert counter == sum( correct.values() + incorrect.values() + unknowns.values() )

	print "depth\ttrue\tfalse\tunknown"
	if scale == 1:
		for l in range( depth + 1 ):
			print "%d\t%d\t%d\t%d" % (l, correct.get( l, 0 ), incorrect.get( l, 0 ), unknowns.get( l, 0 ))
	else:
		for l in range( depth + 1 ):
			print "%d\t%.2f\t%.2f\t%.2f" % (l, scale*correct.get( l, 0 ), scale*incorrect.get( l, 0 ), scale*unknowns.get( l, 0 ))
