#!/bin/bash

globaloptions=$1 #for all plots
shift

options="u $1 w l" #which columns
shift

plotarg="'$1' $options" # first file
shift
files=$*
for file in $files
do
	plotarg=$plotarg", '$file' $options"
done

gnuplot -persist << EOF

#set log 
set nokey
#set yrange [-1:2]

plot $globaloptions $plotarg

EOF
