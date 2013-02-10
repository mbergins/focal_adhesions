plot_server_load_day <- function(data.set,file.name) {
    png(file.name,width=900,height=450,pointsize=24);
    par(bty='n',mar=c(2.8,2.5,0.25,0), mgp=c(1.6,0.5,0),xpd=T,lwd=3)
    plot(data.set,typ='l',ylim=c(0,100),ylab='% Load of Server Max',axes=F,xlab='Time (hours)')
    axis(2);
    axis(1,at=seq(0,1440,by=120),labels=seq(0,24,by=2))
    graphics.off();
}

plot_server_load_week <- function(data.set,file.name) {
    png(file.name,width=900,height=450,pointsize=24);
    par(bty='n',mar=c(2.8,2.5,0.25,0), mgp=c(1.6,0.5,0),xpd=T,lwd=3)
    plot(data.set,typ='l',ylim=c(0,100),ylab='% Load of Server Max',axes=F,xlab='Time (Days)')
    axis(2);
    axis(1,at=seq(0,7*24*60,by=24*60),labels=seq(0,7))
    graphics.off();
}

plot_server_load_all_time <- function(data.set,file.name) {
    png(file.name,width=900,height=450,pointsize=24);
    par(bty='n',mar=c(3,2.5,1.3,0), mgp=c(1.6,0.5,0),xpd=T,lwd=3)
    plot(data.set,typ='l',ylim=c(0,4),ylab='Server Load Level',axes=F,xlab='Time (weeks)')
    axis(2);
    axis(1,at=seq(0,dim(all_time)[1],by=24*60*7),labels=seq(0,dim(all_time)[1]/(24*60*7)))
    graphics.off();
}

args = commandArgs(trailingOnly=T);
if (length(args) != 0) {

	#split out the arguments from the passed in parameters and assign variables 
	#in the current scope
    for (this_arg in commandArgs(TRUE)) {
        split_arg = strsplit(this_arg,"=",fixed=TRUE)
        if (length(split_arg[[1]]) == 1) {
            assign(split_arg[[1]][1], TRUE);
        } else {
            assign(split_arg[[1]][1], split_arg[[1]][2]);
        }
    }
    
    if (! exists('target_dir')) {
        target_dir = '.';
    }

    this_data_set = read.csv(data_file,header=F);
    this_data_set = this_data_set[,1];
    this_data_set[this_data_set > 4] = 4;
	this_data_set = 100*(this_data_set/4);
    switch(plot_type,
           day = plot_server_load_day(this_data_set,file.path(target_dir,'server_load_day.png')),
           week = plot_server_load_week(this_data_set,file.path(target_dir,'server_load_week.png')),
           all_time = plot_server_load_day(this_data_set,file.path(target_dir,'server_load_all_time.png'))
           )
}
