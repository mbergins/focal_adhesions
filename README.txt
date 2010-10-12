Congradulations, you have downloaded the focal adhesion analysis suite. You
should see three folders:

	-data: used to hold the raw imaging data, along with config files
	controlling the options used for processing are stored here
	
	-doc: in here you will find a latex document going into more depth on how
	the series of programs that make up the pipeline are arranged and what each
	program is responsible for

	-src: all the source code that processes the images and conducts the
	statistical analysis is in this folder

In addition to the programs in the src folder, I also expect that you have
working installations of MATLAB and R.

There is a master control script that should be used to process your image,
located at 'src/scripts/build_all_images.pl'. It expects a single command line
parameter, namely a config file. In order to get started a sample data set has
been provided. This zip file should be unzipped into the data directory.
Working from the 'src/scripts' directory, exclude the command:

	./build_all_results.pl -cfg ../../data/config/FA_default.cfg

Now wait for several hours (~5.5 hrs on a 2.67 Ghz 8-core i7 with 12 gigs RAM).
You should see several status messages print to the terminal window. The longest
running part of the script is 'collect_fa_image_set'. To put new data into the
system, I would suggest modeling your data folder after the sample data set
provided, although that can be modified through the config files.

After the programs run, you should see a new directory called results at the top
level of the directory structure. The subfolders within this directory hold all
the image processing results. There are a wide range of temporary files created
during the image processing and analysis steps, some of the output files you
might be interested in are:

	-'visualizations/tracking_seq/edge_track': holds a set of images showing
	each detected adhesion (highlighted in a different color throughout its
	lifetime) and an overlay showing where each of the cell edges were detected
	
	-'visualizations/single_ad': holds sets of small mulitple images showing the
	adhesion of interest in green, other adhesions in blue and the cell edge in
	red, the name of each file corresponds to the lineage number assigned during
	image processing

	-'adhesion_props/assemb_disassem_rates_filtered.csv': This CSV formated file
	contains information about each of the adhesion assembly and disassembly
	rates filtered to only include the significant and positive rates, you can
	see the unfiltered rates in the with a similar name in the same folder. If
	you are interested in seeing a specific adhesion from the experiment, cross
	reference the "lineage_num" with the images in 'visulizations/single_ad'
