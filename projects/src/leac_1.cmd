# set current directory to working directory
#$ -cwd
# combine PBS standard output and error files
#$ -j y

$HOME/git_repos/tracy-3.5-3.2/projects/src/leac param.dat 1
