#!/bin/sh
rm -f rrfileio.h.auto
cat *.h >> rrfileio.h.auto
cp rrfileio.h.auto ../include/rrfileio.h
