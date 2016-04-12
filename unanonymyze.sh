#!/bin/bash

for f in `find . -name '*.m'`; do
	 sed -i 's/org.anon.cocoa/nl.coenvl.sam/g' $f
done
