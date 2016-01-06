#!/bin/sh
FILES="drawcontours.h
drawfield.h
window2d.h"

rm -f rropengl.h.auto
for file in $FILES
do
    if [ ! -e "$file" ]
    then
	echo "$file does not exist."; echo
	continue
    fi
    cat $file >> rropengl.h.auto
done
cp rropengl.h.auto ../include/rropengl.h
