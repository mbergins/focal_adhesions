#!/bin/bash

tmux new-session -d -s FAAS_testing

tmux new-window -t FAAS_testing:0 -n 'testing_server'
tmux split-window -h 

tmux select-pane -t 0
tmux send-keys "tail -f ../logs/development.log" C-m
 
tmux select-pane -t 1
tmux send-keys "cd ../public/" C-m
tmux send-keys "plackup -R ../lib/ ../bin/app.pl" C-m

tmux attach-session -t FAAS_testing
