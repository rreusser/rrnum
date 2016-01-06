#!/bin/sh
rm -f rrgrid.h.auto
cat *.h >> rrgrid.h.auto
cp rrgrid.h.auto ../include/rrgrid.h
