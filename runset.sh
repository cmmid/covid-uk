#!/bin/bash

for i in `seq $2 $3`;
do
    echo "Rscript UK.R $i $1 50"&
done 
