# set current directory to working directory
#$ -cwd
# combine PBS standard output and error files
#$ -j y

$HOME/git_repos/tracy-3.5-3.1/projects/src/ptc/fmap flat_file.dat 1
