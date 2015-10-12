#!/usr/bin/awk -f

/^>/ {
  if( id != "" ) {
    printf "%s\t%s\n", id, sum;
  }
  id=substr( $0, 2 );
  sum = 0;
}

! /^>/ {
  sum+=length($0);
}

END {
  printf "%s\t%s\n", id, sum;
}
