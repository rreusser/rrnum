#!/bin/sh
rm -f rrgeneral.h.auto
cat *.h >> rrgeneral.h.auto
cp rrgeneral.h.auto ../include/rrgeneral.h
