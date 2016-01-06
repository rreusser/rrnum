#!/bin/sh
rm -f rrdata.h.auto
cat *.h >> rrdata.h.auto
cp rrdata.h.auto ../include/rrdata.h
