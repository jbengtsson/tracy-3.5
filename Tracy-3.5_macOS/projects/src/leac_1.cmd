# set current directory to working directory
#$ -cwd
# combine PBS standard output and error files
#$ -j y

$TRACY_LIB/projects/src/leac param.dat 1
