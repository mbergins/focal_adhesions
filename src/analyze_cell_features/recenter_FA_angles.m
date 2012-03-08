function recenter_FA_angles(exp_dir,varargin)

tic;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Setup variables and parse command line
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
i_p = inputParser;

i_p.addRequired('exp_dir',@(x)exist(x,'dir') == 7);
i_p.addParamValue('debug',0,@(x)x == 1 || x == 0);

i_p.parse(exp_dir,varargin{:});

%Add the folder with all the scripts used in this master program
addpath(genpath('..'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Main Program
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

FA_angles = csvread(fullfile(exp_dir,'adhesion_props','lin_time_series','Angle_to_FA_cent.csv'));
FA_cent_pos = csvread(fullfile(exp_dir,'adhesion_props','single_props','Adhesion_centroid.csv'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Find Centroid Direction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cent_pos_diffs = zeros((size(FA_cent_pos,1)-1),2);
for i=2:size(FA_cent_pos,1)
    cent_pos_diffs(i-1,1) = (FA_cent_pos(i,1) - FA_cent_pos(i-1,1));
    cent_pos_diffs(i-1,2) = (FA_cent_pos(i,2) - FA_cent_pos(i-1,2))*-1;
end

mean_cent_movement = mean(cent_pos_diffs);

cent_direction = atan2(mean_cent_movement(2),mean_cent_movement(1))*(180/pi);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Recenter and Save Results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

FA_angles = FA_angles - cent_direction;

FA_angles(FA_angles < -180) = FA_angles(FA_angles < -180) + 360;
FA_angles(FA_angles > 180) = FA_angles(FA_angles > 180) - 360;

csvwrite(fullfile(exp_dir,'adhesion_props','lin_time_series','FA_angle_recentered.csv'),FA_angles);