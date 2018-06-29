#!/bin/sh

prm1=${1-0}

cd $1/dnu
\cp flat_file.fit flat_file.dat
~/git_repos/tracy-3.5_temp/projects/src/get_dnu.sh
~/git_repos/tracy-3.5_temp/projects/src/bare_da.sh
