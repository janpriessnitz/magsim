#!/bin/bash


for simdir in `ls simulations`;
do
    mkdir -p simulations/$simdir/results
    tar -xf simulations/$simdir/results.tar.gz --directory simulations/$simdir/results
done