#!/usr/bin/env bash
echo "Syncing Directories to $1";
rsync --progress --exclude-from .standard_exclude --exclude=**config** -a  $1:~/Documents/Projects/focal_adhesions/data/ ../../data/;
