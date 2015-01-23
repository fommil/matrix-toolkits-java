#!/bin/sh

if [ `git diff | wc -l` -ge 1 ] ; then
    echo "Code formatting does not meet the project's standards:"
    git --no-pager diff
    exit 1
fi
