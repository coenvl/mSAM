#!/bin/bash

DST="../CoCoAnon"

if [ -e "$DST" ]; then
	echo "Target folder $DST already exists"
	exit 1;
fi

mkdir $DST
cp -r functions $DST/.
cp -r scripts $DST/.
cp startup.m $DST/.
cp README.md $DST/.

mkdir $DST/lib
for f in `ls lib/*.jar`; do
    cp $f $DST/`echo $f | sed 's/nl.coenvl.sam/org.anon.cocoa/g'`
done

for f in `find $DST -name '*.m'` $DST/README.md; do
	sed -i 's/nl.coenvl.sam/org.anon.cocoa/g' $f
	sed -i 's/coenvanl@gmail.com/anon187c@gmail.com/g' $f
	sed -i 's/c.j.vanleeuwen-2@tudelft.nl/anon187c@gmail.com/g' $f
	sed -i 's/Coen van Leeuwen/Anomymous/g' $f
	sed -i 's/Coen/Anomymous/g' $f
	sed -i 's/TNO/Anonymous/g' $f
	sed -i 's/SAM/CoCoA/g' $f
done
