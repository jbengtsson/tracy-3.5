# set current directory to working directory
#$ -cwd
# combine PBS standard output and error files
#$ -j y

$HOME/git_repos/projects/src/ptc/fmap flat_file.dat 2
