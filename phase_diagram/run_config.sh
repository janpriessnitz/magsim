#!/bin/bash
d=$1
b=$2
simdir=data/run_${d}_${b}
mkdir -p $simdir
cp metropolis.conf $simdir/
echo "D ${d}e-2" >> $simdir/metropolis.conf
echo "Bz ${b}e-2" >> $simdir/metropolis.conf
echo "#D B" > $simdir/config.info
echo "${d}e-2 ${b}e-2" >> $simdir/config.info
pushd $simdir
../../magsim metropolis.conf
popd
