# this file contains common functions and data for working with confusion
# matrices and their evaluation

from itertools import count
from math import isnan, ceil, log
from sys import stderr

from numpy import array, empty, arange, mean, std, zeros


float_nan = float('NaN')


class ConfusionMatrix:
    def __init__(self, mat, rownames, colnames, title=""):
        self._mat = array(mat)
        self._rownames = rownames
        self._colnames = colnames
        self._rowindex = dict(zip(rownames, count()))
        self._colindex = dict(zip(colnames, count()))
        self.title = title

    # TODO: pre-compute sums to make operations constant in time

    def recall_freqs(self):
        for name, row in zip(self._rownames, self._mat):
            size = row.sum()
            cindex = self._colindex[name]
            try:
                correct = row[cindex]
                if size:
                    yield name, size, correct
            except KeyError:
                pass

    _recall_freqs = recall_freqs

    def recall_freq(self, name):
        try:
            rindex = self._rowindex[name]
            size = self._mat[rindex].sum()
        except KeyError:
            return 0, 0
        try:
            cindex = self._colindex[name]
            correct = self._mat[rindex][cindex]  # TODO: format [a,b]?
        except KeyError:
            correct = 0
        return size, correct

    def _recalls(self):
        for name, size, correct in self._recall_freqs():
            yield name, correct / float(size)

    def recall(self, name=None):
        if name is None:
            return self._recalls()
        size, correct = self.recall_freq(name)
        return correct / float(size)

    def macro_recall(self, ignore_class=""):
        recs = []
        for name, rec in self.recall():
            if name != ignore_class:
                recs.append(rec)
        if recs:
            return mean(recs), std(recs), len(recs)
        return float_nan, float_nan, 0

    # same as accuracy!
    def micro_recall(self, ignore_class=""):
        totalcorrect = 0
        totalsize = 0
        for name, size, correct in self._recall_freqs():
            if name != ignore_class and size:
                totalcorrect += correct
                totalsize += size
        return totalcorrect / float(totalsize)

    def accuracy(self, ignore_class=""):
        totalsize = 0
        totalcorrect = 0
        for name, row in zip(self._rownames, self._mat):
            if name != ignore_class:
                totalsize += row.sum()
                cindex = self._colindex[name]
                try:
                    totalcorrect += row[cindex]
                except KeyError:
                    pass
        if totalsize:
            return totalcorrect / float(totalsize)
        return float_nan

    def entropy_freqs(self, ignore_class=""):
        for i, name in enumerate(self._colnames):
            if name == ignore_class:
                continue
            h_cluster = .0
            col = self._mat[:, i]
            size = float(col.sum())
            for j, rname in zip(col, self._colnames):
                if j and rname != ignore_class:
                    p = j/size
                    h_cluster -= j*log(p, 2.)
            yield name, h_cluster, size  # important: h_cluster is not entropy, but size*entropy

    def entropy(self, ignore_class=""):  # as in doi://10.1093/bioinformatics/btm134
        l = self._mat.shape[0] - (ignore_class in self._rowindex)
        totalsize = .0
        h = .0
        for name, hcl, size in self.entropy_freqs(ignore_class):
            totalsize += size
            # stderr.write("%s has entropy %.2f\n" % (name, hcl))
            h += hcl
        return h/(totalsize*log(l, 2.))

    def misclassification_rate(self, ignore_class=""):
        try:
            ignorecolindex = self._colindex[ignore_class]
        except KeyError:
            ignorecolindex = None
        totalreject = 0
        totalsize = 0
        totalcorrect = 0
        for name, row in zip(self._rownames, self._mat):
            if name != ignore_class:
                totalsize += row.sum()
                cindex = self._colindex[name]
                try:
                    totalcorrect += row[cindex]
                except KeyError:
                    pass
                if ignorecolindex is not None:
                    totalreject += row[ignorecolindex]
        if totalsize:
            return (totalsize - totalcorrect - totalreject) / float(totalsize)
        return float_nan

    def precision_freqs(self):
        for name, i in zip(self._colnames, count()):
            col = self._mat[:, i]
            size = col.sum()
            rindex = self._rowindex[name]
            try:
                correct = col[rindex]
                if size:
                    yield name, size, correct
            except KeyError:
                pass

    _precision_freqs = precision_freqs

    def precision_freq(self, name):
        try:
            cindex = self._colindex[name]
            size = self._mat[:, cindex].sum()
        except KeyError:
            return 0, 0
        try:
            rindex = self._rowindex[name]
            correct = self._mat[rindex][cindex]  # TODO: format [a,b]?
        except KeyError:
            correct = 0
        return size, correct

    def _precisions(self):
        for name, size, correct in self._precision_freqs():
            yield name, correct / float(size)

    def precision(self, name=None):
        if name is None:
            return self._precisions()
        ret = self.precision_freq(name)
        return ret[2] / float(ret[1])

    def macro_precision(self, ignore_class="", truncate=0):
        if not truncate:
            precs = []
            for name, prec in self.precision():
                if name != ignore_class:
                    precs.append(prec)
        else:
            totalsize = 0
            precs_sizes = []
            for name, size, correct in self._precision_freqs():
                if name != ignore_class:
                    totalsize += size
                    precs_sizes.append((size, correct / float(size)))
            precs_sizes.sort(reverse=True)  # reverse sort from high to low bins
            if type(truncate) == float and truncate < 1.:
                precs = []
                lastsize = cumsize = 0
                threshold = ceil(totalsize * truncate)
                for size, prec in precs_sizes:
                    if cumsize > threshold and size < lastsize:  # treat equal size classes
                        break
                    precs.append(prec)
                    cumsize += size
                    lastsize = size
            elif type(truncate) == int:
                precs = [p[1] for p in precs_sizes[:truncate]]
                try:
                    lastsize = precs_sizes[truncate - 1][0]  # TODO: format [a,b]?
                    for size, prec in precs_sizes[truncate - 1:]:
                        if size < lastsize:
                            break
                        precs.append(prec)
                except IndexError:
                    pass
            else:
                raise TypeError(
                    "truncate must either be the number of classes (integer) or a valid fraction (float) between 0 and 1")
        if precs:
            return mean(precs), std(precs), len(precs)
        return float_nan, float_nan, 0

    # same as accuracy!
    def micro_precision(self, ignore_class=""):
        totalcorrect = 0
        totalsize = 0
        for name, size, correct in self._precision_freqs():
            if name != ignore_class and size:
                totalcorrect += correct
                totalsize += size
        return totalcorrect / float(totalsize)
    
    def rand(self, ignore_class=""):  # as in doi://10.1038/nmeth.3103
        npair = lambda i: i*(i-1)  # not including division by 2 because we can often cancel it
        all_pair_sum = 0.
        total_sum = 0.
        row_pair_sum = 0.
        col_sums = zeros(self._mat.shape[1])  # ignore_class entry simply remains zero
        for cname, col in zip(self._rownames, self._mat):
            if cname == ignore_class:
                continue
            row_sum = 0.
            for i, rname, val in zip(count(), self._colnames, col):
                if rname == ignore_class:
                    continue
                all_pair_sum += npair(val)
                row_sum += val
                col_sums[i] += val
            total_sum += row_sum
            row_pair_sum += npair(row_sum)
        col_pair_sum = npair(col_sums).sum()
        # stderr.write("%.2f, %.2f, %.2f\n" % (all_pair_sum, row_pair_sum, col_pair_sum))
        
        total_pair_sum = npair(total_sum)
        if total_pair_sum == 0:
            rand = float('nan')
            arand = float('nan')
        else:
            # simple rand index
            rand = 1. + (2.0*all_pair_sum - row_pair_sum - col_pair_sum)/total_pair_sum
            # adjusted rand index
            t1 = 2*row_pair_sum*col_pair_sum/total_pair_sum
            arand = (all_pair_sum - t1)/((row_pair_sum + col_pair_sum)/2.0 - t1)
        return rand, arand

    def plot_matrix(self, ignore_class="", title="", dpi=300, output=None, fmt=None, extratxt=None, axislabels=True, groupcols=("red", "blue", "grey")):
        import matplotlib
        matplotlib.use('Agg')
        from matplotlib import pyplot, colors, font_manager
        from math import sqrt

        # some constants
        numcolors = 1000
        gridgrey = "0.7"
        gridlwidth = 0.5
        fsize = 10
        labellen = 25
        dampingexp = 1 / 4.
        # labelfont = font_manager.FontProperties( fname="/.../MyriadPro-LightCond.otf" )

        # create temporary color array
        ca = empty((self._mat.shape[0], self._mat.shape[1], 3), dtype=float)  # TODO: not shape[1]?

        # TODO: the right way would be to transform the color map
        denominator = self._mat.sum() ** dampingexp
        transform = lambda l: int(round((l ** dampingexp) * numcolors / denominator))

        # true predictions will be blue, false will be red and rejects will be grey
        reds = colors.LinearSegmentedColormap.from_list("customreds", ["white", groupcols[0]],
                                                        N=numcolors)  # pyplot.get_cmap( "Reds" )
        blues = colors.LinearSegmentedColormap.from_list("customblues", ["white", groupcols[1]],
                                                         N=numcolors)  # pyplot.get_cmap( "Blues" )
        greys = colors.LinearSegmentedColormap.from_list("customgreys", ["white", groupcols[2]],
                                                         N=numcolors)  # pyplot.get_cmap( "Greys" )

        # cut labels for processing
        def cutoff(s):
            if len(s) > labellen:
                return "%s%s" % (s[:labellen], "...")
            return s

        rownames_cut = map(cutoff, self._rownames)
        colnames_cut = map(cutoff, self._colnames)

        # calculate colors for each cell
        for rname, orow, crow in zip(self._rownames, self._mat, ca):
            for cname, i in zip(self._colnames, count()):
                if cname == ignore_class or rname == ignore_class:
                    crow[i] = greys(transform(orow[i]))[:3]
                elif cname == rname:
                    crow[i] = blues(transform(orow[i]))[:3]
                else:
                    crow[i] = reds(transform(orow[i]))[:3]

        # set fonts
        titlefont = font_manager.FontProperties(family="sans-serif", stretch="normal", weight="normal", size="medium",
                                                style="normal")
        extrafont = font_manager.FontProperties(family="serif", stretch="normal", weight="normal", size="small",
                                                style="italic")
        rownum, colnum = self._mat.shape

        labelfont = font_manager.FontProperties(
            family="sans-serif", stretch="condensed", weight="light", size="small", style="normal")

        ax = pyplot.gca()  # fig.add_subplot( 111 )
        im = ax.imshow(ca, interpolation="nearest", origin="upper", aspect="equal")

        # adjust tick fonts
        labelfont_size = (150.0 - min(max(rownum, colnum), 150.0))*0.06  # dynamic adjustment of font size
        # stderr.write("rownum: %i | colum: %i | font size: %.2f\n" % (rownum, colnum, labelfont_size))
        labelfont_y = font_manager.FontProperties(
            family="sans-serif", stretch="condensed", weight="light", size=labelfont_size, style="normal")
        labelfont_x = font_manager.FontProperties(
            family="sans-serif", stretch="condensed", weight="light", size=labelfont_size, style="normal")

        ax.xaxis.set_ticks_position('both')

        ax.xaxis.set_ticks(range(len(self._colnames)))
        ax.xaxis.set_ticks(arange(0.5, len(self._colnames) - 1), minor=True)
        ax.xaxis.set_ticklabels(colnames_cut, fontproperties=labelfont_x)
        # ax.xaxis.grid( True, which="minor", linestyle="-", linewidth=gridlwidth, color=gridgrey )
        ax.xaxis.tick_top()

        ax.yaxis.set_ticks(range(len(self._rownames)))
        ax.yaxis.set_ticks(arange(0.5, len(self._rownames) - 1), minor=True)     
        ax.yaxis.set_ticklabels(rownames_cut, fontproperties=labelfont_y)
        # ax.yaxis.grid( True, which="minor", linestyle="-", linewidth=gridlwidth, color=gridgrey )
        ax.yaxis.tick_left()

        for label in ax.xaxis.get_ticklabels():
            label.set_rotation(90)
            label.set_fontproperties(labelfont_x)
        
        if not axislabels or labelfont_size < 0.5:  # minimal font size threshold
            ax.set_xticklabels('', visible=False)
            ax.set_yticklabels('', visible=False)

        if title:
            # ax.set_title( title )
            ax.set_xlabel(title, labelpad=10, fontproperties=titlefont)
        else:
            ax.set_xlabel(self.title, labelpad=10, fontproperties=titlefont)
        # ax.set_title( self._title )

        ax.tick_params(which="minor", direction="out", length=labelfont_size, width=gridlwidth, color=gridgrey)
        ax.tick_params(which="major", length=0)
        # ax.title.set_y(1.05)
        # ax.margins( 0.1 )

        # shrink height of plot to free space for extratxt
        box = ax.get_position()
        ax.set_position([box.x0, box.y0 + box.height * 0.1, box.width, box.height * 0.9])

        fig = pyplot.gcf()
        # fig.set_size_inches( 10, 10 )
        fig.text(0.1, 0, extratxt, fontproperties=extrafont)
        # fig.tight_layout()

        # pyplot.draw()
        # fig.subplots_adjust( left=1, right=1, top=1, bottom=1 )
        # pyplot.sci( im ) #probably not needed
        # ax.annotate(extratxt, xy=(0, 0), xycoords="axes points")
        # pyplot.colorbar()
        # pyplot.tight_layout()

        if output:
            pyplot.savefig(output, format=fmt, dpi=dpi, transparent=True, bbox_inches="tight", pad_inches=1.3)
        else:
            pyplot.show()
        pyplot.close()

    plotMatrix = plot_matrix  # TODO: consistent naming scheme (underscore); transitional for backwards compatibility

    def items(self):
        for rname in self._rownames:
            for cname in self._colnames:
                yield (rname, cname), self._mat[self._rowindex[rname], self._colindex[cname]]

    def get_col(self, key):  # untested; necessary?
        try:
            return [i for i in self._mat[:, self._colindex[key]]]
        except KeyError:
            return (0,) * len(self._colindex)

    def get_cols(self):  # untested; necessary?
        for cname in self._colnames:
            yield cname, self.get_col(cname)

    def get_row(self, key):  # untested; necessary?
        try:
            return [i for i in self._mat[self._rowindex[key]]]
        except KeyError:
            return (0,) * len(self._rowindex)

    def get_rows(self):  # untested; necessary?
        for rname in self._rownames:
            yield rname, self.get_row(rname)


def parse_confusion_matrix(lines):
    line = None
    try:
        while True:
            # eat white space and comments
            while not line or line == "\n" or line[0] == '#':
                line = lines.next()

            line = line.rstrip("\n").split("\t")
            title, colnames = line[0], line[1:]

            rows = []
            rownames = []

            for line in lines:
                if line == "\n":  # stop at empty line
                    break
                line = line.rstrip("\n").split("\t")
                rownames.append(line[0])
                rows.append(map(float, line[1:]))

            # construct numpy matrix
            # cmat = array( rows )
            yield ConfusionMatrix(rows, rownames, colnames, title)
            line = lines.next()

    except StopIteration:
        pass  # raises StopIteration automatically

# TODO: consistent naming scheme (underscore); transitional for backwards compatibility
parseConfusionMatrix = parse_confusion_matrix
