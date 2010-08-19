function turnover(depth, time_interval)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% the inputs are 
 
%  inputs:        depth: this is the depth of the windows in pixel
%                 time_iterval: this is the time interval for sampling
%
% outputs:        It plots the time variation of 
%                 average poly, depoly, and average turnover
%                 (poly+depoly) events as well as the corresponding normalized 
%                 average values per unit area of the region of interest 
%                 (lamellum or Lamellipodium) in units of (Micron)^2.
%                 
%         
%
% by: Mohsen Sabouri Ghomi 6/12/07
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pixel=67; 

load('segment_length_av');
load ('score');


num_segment=size(score,1);
segment_depth=depth*pixel/1000;
area=segment_length_av*num_segment*segment_depth;
turnover=score;
turnover_per_area=score/area;

poly_ind=turnover>0;
poly=poly_ind.*turnover;
depoly_ind=turnover<0;
depoly=depoly_ind.*turnover;
av_poly=mean(poly);
av_depoly=mean(depoly);
av_turnover=mean(turnover);

poly_per_area_ind=turnover_per_area>0;
poly_per_area=poly_per_area_ind.*turnover_per_area;
depoly_per_area_ind=turnover_per_area<0;
depoly_per_area=depoly_per_area_ind.*turnover_per_area;
av_poly_per_area=mean(poly_per_area);
av_depoly_per_area=mean(depoly_per_area);
av_turnover_per_area=mean(turnover_per_area);


score_x_time_axis = (1:size(score,2))*time_interval;

std_depoly=std(depoly,0,1);
std_poly=std(poly,0,1);
std_turnover=std(turnover,0,1);

std_depoly_per_area=std(depoly_per_area,0,1);
std_poly_per_area=std(poly_per_area,0,1);
std_turnover_per_area=std(turnover_per_area,0,1);

% exist_dir = exist('normalized_turnover_figures', 'dir');
% if exist_dir == 1;
% rmdir normalized_turnover_figures;
mkdir('figures_turnover');
figures_dir = ['figures_turnover' filesep];

%%%%%%%%%%%%%%%%%%  Unnormalized turnover plots   %%%%%%%%%%%%%%%%%%%%%%%%%

poly_unnormalized = figure;
plot(score_x_time_axis, av_poly,'-ro');
hold on;
plot(score_x_time_axis, av_poly+std_poly,'--');
plot(score_x_time_axis, av_poly-std_poly,'--');
h = legend('data','std','std', 'Location','Best');
set(h,'Interpreter','none')
xlabel('Time (s) ');
ylabel('Average poly ( score )  ');
hgsave(poly_unnormalized, [figures_dir, 'av_poly.fig']);

depoly_unnormalized = figure;
plot(score_x_time_axis, av_depoly,'-ro');
hold on;
plot(score_x_time_axis, av_depoly+std_depoly,'--');
plot(score_x_time_axis, av_depoly-std_depoly,'--');
h = legend('data','std','std', 'Location','Best');
set(h,'Interpreter','none')
xlabel('Time (s)');
ylabel('Average depoly ( score ) ');
hgsave(depoly_unnormalized,[figures_dir, 'av_depoly.fig']);

turnover_unnormalized = figure;
plot(score_x_time_axis, av_turnover,'-ro');
hold on;
plot(score_x_time_axis, av_turnover+std_turnover,'--');
plot(score_x_time_axis, av_turnover-std_turnover,'--');
h = legend('data','std','std', 'Location','Best');
set(h,'Interpreter','none')
xlabel('Time (s) ');
ylabel('Average turnover ( score )  '); 
hgsave(turnover_unnormalized, [figures_dir, 'av_turnover.fig']);

%%%%%%%%%%%%%%%%%%%%%   normalized turnover plots   %%%%%%%%%%%%%%%%%%%%%%%

poly_normalized = figure;
plot(score_x_time_axis, av_poly_per_area,'-ro');
hold on;
plot(score_x_time_axis, av_poly_per_area+std_poly_per_area,'--');
plot(score_x_time_axis, av_poly_per_area-std_poly_per_area,'--');
h = legend('data','std','std', 'Location','Best');
set(h,'Interpreter','none')
xlabel('Time (s) ');
ylabel('Average poly (score / \mu m^2)  ');
hgsave(poly_normalized, [figures_dir, 'av_poly_per_area.fig']);

depoly_normalized = figure;
plot(score_x_time_axis, av_depoly_per_area,'-ro');
hold on;
plot(score_x_time_axis, av_depoly_per_area+std_depoly_per_area,'--');
plot(score_x_time_axis, av_depoly_per_area-std_depoly_per_area,'--');
h = legend('data','std','std', 'Location','Best');
set(h,'Interpreter','none')
xlabel('Time (s)');
ylabel('Average depoly (score / \mu m^2) ');
hgsave(depoly_normalized,[figures_dir, 'av_depoly_per_area.fig']);

turnover_normalized = figure;
plot(score_x_time_axis, av_turnover_per_area,'-ro');
hold on;
plot(score_x_time_axis, av_turnover_per_area+std_turnover_per_area,'--');
plot(score_x_time_axis, av_turnover_per_area-std_turnover_per_area,'--');
h = legend('data','std','std', 'Location','Best');
set(h,'Interpreter','none')
xlabel('Time (s) ');
ylabel('Average turnover (score / \mu m^2)  '); 
hgsave(turnover_normalized, [figures_dir, 'av_turnover_per_area.fig']);

end % end turnover 
