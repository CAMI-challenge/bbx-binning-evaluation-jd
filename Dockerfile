FROM bioboxes/base
MAINTAINER Johannes Dröge, johannes.droege@uni-duesseldorf.de

# add required Debian packages here (uncomment)
RUN ${BBX_BINDIR}/dockerfile-install-packages python2.7-numpy python-matplotlib

# add task definitions
COPY tasks ${BBX_TASKDIR}

# add required program files
COPY opt ${BBX_OPTDIR}
