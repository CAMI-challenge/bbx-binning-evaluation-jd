FROM bioboxes/base
MAINTAINER Johannes Dröge, johannes.droege@uni-duesseldorf.de

# add required Debian packages here (uncomment)
#RUN ${BBX_BINDIR}/dockerfile-install-packages mypackage1 mypackage2

# add task definitions
COPY tasks ${BBX_TASKDIR}

# add required program files
COPY opt ${BBX_OPTDIR}

