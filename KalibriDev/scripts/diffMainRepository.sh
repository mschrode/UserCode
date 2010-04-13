#!/bin/bash

for file in `ls -1 *.h *.cc`; do
#    echo $file
    diff --brief -I '^//!*\s*\WId: .*\$' $file ~/Kalibri/$file > /dev/null 2>&1
    if [ "$?" == 1 ]; then
	if [ "$1" == "verbose" ]; then
	    echo ""
	    echo ""
	fi
	echo "Differences in file $file"
	if [ "$1" == "verbose" ]; then
	    diff -I '^//!*\s*\WId: .*\$' $file ~/Kalibri/$file
	fi
    elif [ "$?" == 2 ]; then
	echo "File $file does not exist in ~/Kalibri"
    fi
done
