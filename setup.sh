#!/bin/sh

# Check for dependencies
if ! hash git 2>/dev/null; then
    echo "git must be installed to run automatic setup!"
fi
if ! hash curl 2>/dev/null; then
    echo "wget must be installed to run automatic setup!"
fi
if ! hash matlab 2>/dev/null; then
    echo "Matlab must be installed and be available in PATH (as 'matlab') to run notebooks!"
fi

# Change current directory to the directory of this script
dir=$(cd -P -- "$(dirname -- "$0")" && pwd -P)
cd $dir
julia -e 'using Pkg; Pkg.instantiate()'
# change to parent directory
cd ../
# clone required repositories
git clone https://github.com/JeffFessler/mirt.git
git clone https://github.com/hakkelt/reproduce-l-s-dynamic-mri.git
# download data
cd reproduce-l-s-dynamic-mri
mkdir data
cd data
curl -O https://web.eecs.umich.edu/~fessler/irt/reproduce/19/lin-19-edp/data/Xinf.mat
curl -O https://web.eecs.umich.edu/~fessler/irt/reproduce/19/lin-19-edp/data/abdomen_dce_ga.mat
curl -O https://web.eecs.umich.edu/~fessler/irt/reproduce/19/lin-19-edp/data/aperiodic_pincat.mat
curl -O https://web.eecs.umich.edu/~fessler/irt/reproduce/19/lin-19-edp/data/cardiac_cine_R6.mat
curl -O https://web.eecs.umich.edu/~fessler/irt/reproduce/19/lin-19-edp/data/cardiac_perf_R8.mat
curl -O https://web.eecs.umich.edu/~fessler/irt/reproduce/19/lin-19-edp/data/readme-data.txt

echo "Automatic setup completed!"