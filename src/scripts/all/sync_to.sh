#!/usr/bin/env bash
echo "Syncing Directories to $1";
rsync --exclude-from .standard_exclude -a ../../ $1:~/Documents/Projects/focal_adhesions/;
