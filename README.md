# A biobox for CAMI binning evaluation

## Simple instructions
      docker build -t="myuser/binning-evaluation" .
      cd "$(mktemp -d)"
      mkdir input output cache
      # now place biobox.yaml and other files into input folder
      docker run -v $PWD/input:/bbx/mnt/input:ro -v $PWD/output:/bbx/mnt/output:rw -v $PWD/cache:/bbx/mnt/cache:rw myuser/binning-evaluation

