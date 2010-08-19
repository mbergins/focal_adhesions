#!/usr/bin/env bash
echo "Syncing Config to $1";
rsync --progress -a ../../data/config/ $1:~/Documents/Projects/focal_adhesions/data/config;
