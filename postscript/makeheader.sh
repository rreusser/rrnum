#!/bin/sh
FILES="drawpscontours.h"

rm -f rrps.h.auto
for file in $FILES
do
    if [ ! -e "$file" ]
    then
	echo "$file does not exist."; echo
	continue
    fi
    cat $file >> rrps.h.auto
done
cp rrps.h.auto ../include/rrps.h
