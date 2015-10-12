#!/usr/bin/env python
# This script takes two files, the first containing class labels and the second
# containing class predictions in the following format:
#
# entity_name\tclass_name_A\tclass_name_B...
#
# where \t means TAB separation and after the entity name there can be as many
# columns as you which but at least one.
#
# For each column beginning from the second column, the script will calculate
# the confusion matrix (contingency table) containing the frequencies of
# combinations between label and predicted classes and write them to standard
# output.
#
# Conventions:
# a) Strict one-to-one line order is not neccessary
# b) If predictions are missing, script will exit and tell you (override params).
# c) If labels are missing, script will exit and tell you.
# d) Output order of rows and columns is alphabetical
# e) Comment lines in input must start with "#" (first character)

from sys import argv, stdout, stderr, stdin, exit
from itertools import repeat, product

# TODO: add missing predictions to reject class (ignore_class -> reject_class)
# TODO: limit to classes constituing at least x % of sample sequences

# print information on usage
def usage():
    print >> stderr, "Usage: ", argv[
        0], "--rows label.racol --columns predictions.racol [ --weights seq.length --matrix-form sparse/quadratic" \
            " --class-for-missing-predictions "" --allow-missing-rows --allow-missing-columns" \
            " --multiclass-separator ';']"

# helper function
def quadratic_axes(s1, s2, typeconv=None):
    classes = list(s1 | s2)
    if typeconv:
        classes = map(typeconv, classes)
        classes.sort()
        classes = map(str, classes)
    else:
        classes.sort()
    return classes, classes

# helper function
def sparse_axes(s1, s2, typeconv=None):
    classes1 = list(s1)
    classes2 = list(s2)
    if typeconv:
        classes1 = map(typeconv, classes1)
        classes2 = map(typeconv, classes2)
        classes1.sort()
        classes2.sort()
        classes1 = map(str, classes1)
        classes2 = map(str, classes2)
    else:
        classes1.sort()
        classes2.sort()
    return classes1, classes2

# helper function
swapzip = lambda a, b: zip(b, a)

# simple dummy weight function counting each sequences as one
class oneweight:
    __getitem__ = lambda self, key: 1


if __name__ == "__main__":
    import getopt

    # parse command line options
    try:
        opts, args = getopt.getopt(argv[1:], "h1:2:w:t:m:c:bans:",
                                   ["help", "rows=", "columns=", "weights=", "title=", "matrix-form=",
                                    "class-for-missing-columns=", "allow-missing-rows",
                                    "allow-missing-columns", "numeric-classes", "multiclass-separator="])
    except getopt.GetoptError, err:
        # print help information and exit:
        print str(err) # will print something like "option -a not recognized"
        usage()
        exit(2)

    # defaults
    matrix_form = "sparse"
    axes = sparse_axes
    label_filename = None
    pred_filename = "-"
    weight_filename = None
    relax_matches = False
    class_for_missing_predictions = ""
    allow_missing_labels = False
    allow_missing_predictions = False
    title = ""
    multiclass_separator = ""
    numeric_classes = False

    # option parsing
    for o, a in opts:
        if o in ("-h", "--help"):
            usage()
            exit()
        elif o in ("-1", "--rows"):
            label_filename = a
        elif o in ("-2", "--columns"):
            pred_filename = a
        elif o in ("-w", "--weights"):
            weight_filename = a
        elif o in ("-t", "--title"):
            title = a
        elif o in ("-m", "--matrix-form"):
            if a == "quadratic":
                axes = quadratic_axes
                matrix_form = "quadratic"
            elif a == "sparse":
                print >> stderr, "setting sparse matrix"
                axes = sparse_axes
                matrix_form = "sparse"
        elif o in ("-c", "--class-for-missing-predictions"):
            class_for_missing_predictions = a
        elif o in ("-b", "--allow-missing-rows"):
            allow_missing_labels = True
        elif o in ("-a", "--allow-missing-columns"):
            allow_missing_predictions = True
        elif o in ("-s", "--multiclass-separator"):
            multiclass_separator = a
        elif o in ("-n", "--numeric-classes"):
            numeric_classes = True
        else:
            assert False, "unhandled option"

    # check parameters
    if not label_filename:
        print >> stderr, "you must specify a label file"
        usage()
        exit(3)

    # be verbose
    print >> stderr, "Using row file", label_filename
    print >> stderr, "Using column file", pred_filename
    print >> stderr, "Using matrix format", matrix_form
    print >> stderr, "Using multiclass separator", multiclass_separator

    # read weights if given
    if weight_filename:
        print >> stderr, "Using weight file", weight_filename
        weight = {}
        with open(weight_filename, "r") as f:
            for line in f:
                name, w = line.strip().split("\t", 2)[:2]
                weight[name] = float(w)
    else:
        weight = oneweight()

    # before doing the actual processing we peek into the label file to find out
    # the number of columns
    firstline = open(label_filename, "r").next()
    num_tables = firstline.count("\t")
    if firstline[0] == "#":
        titles = firstline[1:].strip().split("\t")[1:]
    elif title:
        titles = [title for i in xrange(num_tables)]
    else:
        titles = ["confusion matrix" for i in xrange(num_tables)]
    tables = [{} for i in xrange(num_tables)]
    class_names_per_table = [( set(), set() ) for i in xrange(num_tables)]

    print >> stderr, "Calculating %i confusion matrices from input" % (num_tables)

    label_filehandle = open(label_filename, "r")

    if pred_filename == "-":
        pred_filehandle = stdin
    else:
        pred_filehandle = open(pred_filename, "r")

    # data holder objects
    cache = ({}, {})
    not_empty = [True, True]

    while any(not_empty):
        process_pairs = [] #local workload
        process_weights = []
        for fhandle, index_this, index_other, zipfnct in zip((label_filehandle, pred_filehandle), (0, 1), (1, 0),
                                                             (zip, swapzip)):
            if not_empty[index_this]:
                try:
                    line = fhandle.next()
                    if line[0] != "#":
                        line = line.rstrip("\n").split("\t")
                        name, classes = line[0], line[1:]

                        try:
                            classes_cached = cache[index_other].pop(name) #look in cache
                            process_pairs.append(zipfnct(classes, classes_cached))
                            process_weights.append(weight[name])
                        except KeyError:
                            cache[index_this][name] = classes #put into cache

                except StopIteration:
                    not_empty[index_this] = False

        # sum up frequencies in tables
        for pairs_list, w in zip(process_pairs, process_weights):
            for table, pair, class_names in zip(tables, pairs_list, class_names_per_table):
                if multiclass_separator:
                    pair = product(pair[0].split(multiclass_separator), pair[1].split(multiclass_separator))
                else:
                    pair = (pair,)

                for single_pair in pair:
                    try:
                        table[single_pair] += w
                    except KeyError:
                        table[single_pair] = w

                    # maintain set of all classes for rows/columns each
                    class_names[0].add(single_pair[0])
                    class_names[1].add(single_pair[1])

    label_filehandle.close()
    pred_filehandle.close()

    # check for correct matches
    print >> stderr, len(cache[1]), "labeled entries are missing in label file"
    if not allow_missing_labels and cache[1]:
        print >> stderr, "Not allowed!"
        exit(5)

    print >> stderr, len(cache[0]), "predicted entries are missing in prediction file"
    if cache[0]:
        if not allow_missing_predictions:
            print >> stderr, "Not allowed!"
            exit(6)
        else:
            # adding missing predictions to class_for_missing_predictions
            c = 0
            for name, classes in cache[0].items():
                w = weight[name]
                for table, pair, class_names in zip(tables, zip(classes, repeat(class_for_missing_predictions)),
                                                    class_names_per_table):
                    # support multiple classes with separator
                    if multiclass_separator:
                        pair = (pair[0].split(multiclass_separator), pair[1].split(multiclass_separator))
                    for single_pair in product(*pair):
                        try:
                            table[single_pair] += w
                        except KeyError:
                            table[single_pair] = w

                        # maintain set of all classes for rows/columns each
                        class_names[0].add(single_pair[0])
                        class_names[1].add(single_pair[1])
                c += 1
            print >> stderr, "%i missing predictions were put into class with name \"%s\"" % (
            c, class_for_missing_predictions)

    # print each of the tables
    while tables:
        table = tables.pop()
        classes = class_names_per_table.pop()
        title = titles.pop()
        if numeric_classes:
            row_classes, column_classes = axes(*classes, typeconv=int)
        else:
            row_classes, column_classes = axes(*classes)

        # print row header
        stdout.write("%s\t" % (title))
        print "\t".join(column_classes)

        for crow in row_classes:
            #print row header entry
            stdout.write("%s\t" % (crow))

            # print entries in defined order
            for ccol in column_classes[:-1]:
                stdout.write("%.2f\t" % (table.get((crow, ccol), 0)))
            stdout.write("%.2f\n" % (table.get((crow, column_classes[-1]), 0)))
        print "\n"
