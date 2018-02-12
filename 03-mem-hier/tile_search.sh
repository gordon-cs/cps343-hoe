#!/bin/bash
#----------------------------------------------------------------------------
#----------------------------------------------------------------------------
# Generate data for matrix-matrix product C = AB
#----------------------------------------------------------------------------
#----------------------------------------------------------------------------

#----------------------------------------------------------------------------
# Function to multiply two floating point values and round result to the
# nearest integer
function mult_and_round()
{
    echo $1 $2 | awk '{print int($1*$2+0.5)}'
}

#----------------------------------------------------------------------------

# Process command line arguments
if [ $# != 3 ] && [ $# != 4 ]
then
    echo "usage: $0 MODE NROWS NCOLS [FACTOR]"
    echo
    echo "where"
    echo -e "\tMODE    is one of ijk, ikj, jik, jki, kij, kji"
    echo -e "\tNROWS   is the number of rows and columns of the result matrix"
    echo -e "\tNCOLS   is the number of columns in the first multiplicand"
    echo -e "\tFACTOR  is the tile size growth factor: (default is 1.3)"
    exit
fi
m=$1        # first argument is string of characters i, j, and k in any order
M=$2        # second argument is number of rows and columns in square result C
N=$3        # third argument is number of columns in A and number of rows in B
g=$4        # optional 4th argument is tile dimension growth factor
g=${g:-1.3} # default if g is empty is 1.3

# loop over exponentially growing tile sizes.  Results are output in two
# columns:
#   tile_size    time_in_seconds
#
for ((ts=36; ts<=M*M; ts=$(mult_and_round $g $ts)))
do
    echo "$ts $(./matrix_prod -r $m $M $N $M $ts | awk '{print $4}')"
done
