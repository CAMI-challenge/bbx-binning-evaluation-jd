# cami-binning-eval

A biobox for CAMI binning evaluation

## Simple instructions
      docker build -t="myuser/bbx-eval-binning" .
      cd "$(mktemp -d)"
      mkdir input output cache
      # now place biobox.yaml and other files into input folder
      docker run -v $PWD/input:/bbx/mnt/input:ro -v $PWD/output:/bbx/mnt/output:rw -v $PWD/cache:/bbx/mnt/cache:rw fungs/bbx-eval-binning-joh

