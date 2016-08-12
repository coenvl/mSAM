#!/bin/bash

DST="../CoCoAnon"
GIT="anonymousgit:anonc187/mCoCoA.git"

if [ -e "$DST" ]; then
	echo "Target folder $DST already exists"
	exit 1;
fi

# Clone git repository
git clone $GIT $DST
git -C $DST config user.name "Anon Y. Mous"
git -C $DST config user.email anon187c@gmail.com

echo "Cleaning up old repository"
rm -rf $DST/*

echo "Copying working folder to $DST"
# Copy current working files
cp -r functions $DST/.
cp -r scripts $DST/.
cp startup.m $DST/.
cp README.md $DST/.

mkdir $DST/lib
for f in `ls lib/*.jar`; do
    cp $f $DST/`echo $f | sed 's/nl.coenvl.sam/org.anon.cocoa/g'`
done

echo "Anonymizing $DST"
for f in `find $DST -name '*.m'` $DST/README.md; do
	sed -i 's/nl.coenvl.sam/org.anon.cocoa/g' $f
	sed -i 's/coenvanl@gmail.com/anon187c@gmail.com/g' $f
	sed -i 's/c.j.vanleeuwen-2@tudelft.nl/anon187c@gmail.com/g' $f
	sed -i 's/Coen van Leeuwen/Anomymous/g' $f
    sed -i 's/leeuwencjv/Anomymous/g' $f
	sed -i 's/Coen/Anomymous/g' $f
	sed -i 's/TNO/Anonymous/g' $f
	sed -i 's/SAM/CoCoA/g' $f
done

echo "Sending commits"
cd $DST
git add *
git commit -m `date +%Y-%m-%d`
git push origin master
