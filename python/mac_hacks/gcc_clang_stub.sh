#!/bin/sh

CC=g++-11

COLLECTED_ARGS=""
for i in $*
do
    case ${i} in
	-stdlib=libc++)
	# ignore argument
	;;
	*)
	    COLLECTED_ARGS="$COLLECTED_ARGS $i"
	;;
    esac
done

echo $CC $COLLECTED_ARGS
