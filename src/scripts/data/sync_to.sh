#!/usr/bin/env bash
echo "Syncing Directories to $1";
rsync --progress --exclude-from .standard_exclude --exclude=**config** -a ../../data/ $1:~/Documents/Projects/focal_adhesions/data/;
