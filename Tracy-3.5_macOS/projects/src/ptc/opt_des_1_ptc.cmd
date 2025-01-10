# set current directory to working directory
#$ -cwd
# combine PBS standard output and error files
#$ -j y

$HOME/git_repos/tracy-3.5_temp/projects/src/ptc/opt_des_1_ptc $lat_file
