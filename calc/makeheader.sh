#!/bin/sh
FILES="ab2.h
euler.h
nder.h
rk2.h
stencil.h
integrate.h
nr.h
rk4.h"

rm -f rrcalc.h.auto
for file in $FILES
do
    if [ ! -e "$file" ]
    then
	echo "$file does not exist."; echo
	continue
    fi
    cat $file >> rrcalc.h.auto
done
cp rrcalc.h.auto ../include/rrcalc.h
