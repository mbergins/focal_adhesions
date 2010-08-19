#!/usr/bin/env bash
echo "Syncing Focal Adhesion Results to $1";
rsync --progress -a ../../results/focal_adhesions/ $1:~/Documents/Projects/focal_adhesions/results/focal_adhesions/;

echo "Syncing Sampled Focal Adhesion Results to $1";
rsync --progress -a ../../results/focal_adhesions_reduced/ $1:~/Documents/Projects/focal_adhesions/results/focal_adhesions_reduced/;

echo "Syncing S178A Results to $1";
rsync --progress -a ../../results/S178A/ $1:~/Documents/Projects/focal_adhesions/results/S178A;

echo "Syncing Sampled S178A Results to $1";
rsync --progress -a ../../results/S178A_reduced/ $1:~/Documents/Projects/focal_adhesions/results/S178A_reduced;

echo "Syncing Rap-Src Results to $1";
rsync --progress -a ../../results/rap_src/ $1:~/Documents/Projects/focal_adhesions/results/rap_src/;

echo "Syncing Rap-Src Control Results to $1";
rsync --progress -a ../../results/rap_src_control/ $1:~/Documents/Projects/focal_adhesions/results/rap_src_control/;

echo "Syncing Simulation Results to $1";
rsync --progress -a ../../results/simulation/ $1:~/Documents/Projects/focal_adhesions/results/simulation/;
