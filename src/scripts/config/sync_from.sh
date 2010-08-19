#!/usr/bin/env bash
echo "Syncing Config from $1";
rsync --progress -a $1:~/Documents/Projects/focal_adhesions/data/config/ ../../data/config;
