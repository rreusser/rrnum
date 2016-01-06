#!/bin/sh
FILES="sparse.h
bicg.h
tdiag.h"

rm -f rrlinalg.h.auto
for file in $FILES
do
    if [ ! -e "$file" ]
    then
	echo "$file does not exist."; echo
	continue
    fi
    cat $file >> rrlinalg.h.auto
done
cp rrlinalg.h.auto ../include/rrlinalg.h
