#!/usr/bin/perl -w

system("scp -r mimir2:/var/www/focal_adhesions/data/$ARGV[0] ../../data/");

system("scp -r mimir2:/var/www/focal_adhesions/src/webapp/public/results/$ARGV[0].zip ../../results/");

system("unzip -q ../../results/$ARGV[0].zip -d ../../results/");
system("rm ../../results/$ARGV[0].zip");
