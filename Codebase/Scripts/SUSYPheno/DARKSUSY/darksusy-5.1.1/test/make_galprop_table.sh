#!/bin/sh

# This should usually be invoked as:
#
# sh make_galprop_table.sh ../share/DarkSUSY/FITS/GF_50p_RACC
#   for option 'default'
#
#   or
#
# sh make_galprop_table.sh ../share/DarkSUSY/FITS/GF_50p_DIFF
#   for option 'plain'

for i in {1..51}; do
    ./dstest-galprop-one <<EOF
$i
EOF
done
cat ${1}0*.dat > ${1}foo.dat
