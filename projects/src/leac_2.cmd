# set current directory to working directory
#$ -cwd
# combine PBS standard output and error files
#$ -j y

$TRACY_LIB/tracy-3.5_temp/projects/src/leac param.dat 2
