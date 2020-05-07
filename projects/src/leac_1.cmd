# set current directory to working directory
#$ -cwd
# combine PBS standard output and error files
#$ -j y

$TRACY_DIR/projects/src/leac param.dat 1
