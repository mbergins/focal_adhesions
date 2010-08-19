function ratio_activity(edge_parameters, merg_parameters, post_parameters, cell_variables)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%           Calculate the ratio between activity 1 and 2                  %
%                                                                         %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PR_ALPHA_DIR        = post_parameters.pr_alpha_dir;  
PROJECT_DIR         = post_parameters.project_dir;  
POST_PROC_OPTIONS   = post_parameters.post_proc_options;  
DO_PROT             = post_parameters.do_prot;  
DO_RETRO            = post_parameters.do_retro;  
DO_SCORE            = post_parameters.do_score;  
DO_ACTIVITY_1       = post_parameters.do_activity_1;  
DO_ACTIVITY_2       = post_parameters.do_activity_1;  
DO_DETREND          = post_parameters.do_detrend;
TREND_METHOD        = post_parameters.trend_method;
trend               = post_parameters.trend;
windowSize          = post_parameters.windowSize;
spline_s            = post_parameters.spline_s;
RESULT_DIR          = post_parameters.result_dir;


if DO_PROT
    protrusion = cell_variables.protrusion;
end
if DO_SCORE
    score = cell_variables.score;
end
if DO_RETRO
    retrograde_flow = cell_variables.retrograde_flow;
    shift_retrograde_flow = cell_variables.shift_retrograde_flow;
end
if DO_ACTIVITY_1 | DO_ACTIVITY_2
   activity = cell_variables.activity;
   shift_activity = cell_variables.shift_activity;
end

if DO_PROT
    %get the number of segments
    n_segments = size(protrusion, 1);
    %get the number of time steps
    n_time     = size(protrusion, 2);
else
    %get the number of segments
    n_segments = size(activity, 2);
    %get the number of time steps
    n_time     = size(activity, 3);
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Detrend data               %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if TREND_METHOD > 1
    if DO_PROT
        [protrusion, protrusion_trend, av_de_trend_protrusion_time] = analyzeTimeSeries(protrusion,TREND_METHOD, windowSize, spline_s);
    end
    if DO_ACTIVITY_1
        [activity_1, activity_trend_1, av_de_trend_activity_time_1] = analyzeTimeSeries(squeeze(activity(1,:,:)), TREND_METHOD, windowSize, spline_s);
        [shift_activity_1, shift_activity_trend_1, av_de_trend_shift_activity_time_1] = analyzeTimeSeries(squeeze(shift_activity(1,:,:)),TREND_METHOD, windowSize, spline_s);
    end
    if DO_ACTIVITY_2
        [activity_2, activity_trend_2, av_de_trend_activity_time_2] = analyzeTimeSeries(squeeze(activity(2,:,:)), TREND_METHOD, windowSize, spline_s);
        [shift_activity_2, shift_activity_trend_2, av_de_trend_shift_activity_time_2] = analyzeTimeSeries(squeeze(shift_activity(2,:,:)),TREND_METHOD, windowSize, spline_s);
    end    
    if DO_RETRO
        [retrograde_flow, retrograde_flow_trend, av_de_trend_retrograde_flow_time] = analyzeTimeSeries(retrograde_flow,TREND_METHOD, windowSize, spline_s);
        [shift_retrograde_flow, shift_retrograde_flow_trend, av_de_trend_shift_retrograde_flow_time] = analyzeTimeSeries(shift_retrograde_flow,TREND_METHOD, windowSize, spline_s);
    end
    if DO_SCORE
        [score, score_trend, av_de_trend_score_time] = analyzeTimeSeries(score,TREND_METHOD, windowSize, spline_s);
    end

    if trend == 1
        if DO_PROT
            protrusion = protrusion_trend;
        end
        if DO_ACTIVITY_1
            activity_1 = activity_trend_1;
            shift_activity_1 = shift_activity_trend_1;
        end
        if DO_ACTIVITY_2
            activity_2 = activity_trend_2;
            shift_activity_2 = shift_activity_trend_2;
        end        
        if DO_RETRO
            retrograde_flow = retrograde_flow_trend;
            shift_retrograde_flow = shift_retrograde_flow_trend;
        end
        if DO_SCORE
            score = score_trend;
        end
    end
else
    if DO_ACTIVITY_1
        activity_1 = squeeze(activity(1,:,:));
        shift_activity_1 = squeeze(shift_activity(1,:,:));
        % check if there is only one time step
        if n_time == 1
            activity_1 = activity_1';
            shift_activity_1 = shift_activity_1';
        end
    end  
    if DO_ACTIVITY_2
        activity_2 = squeeze(activity(2,:,:));
        shift_activity_2 = squeeze(shift_activity(2,:,:));
        % check if there is only one time step
        if n_time == 1
            activity_2 = activity_2';
            shift_activity_2 = shift_activity_2';
        end
    end     
end


if DO_ACTIVITY_1 & DO_ACTIVITY_2
    ratio       = activity_1./activity_2;
    shift_ratio = shift_activity_1./shift_activity_2;
    
    % segment average ratio 
    ratio_seg_av = sum(ratio,1) ./ n_segments;
    shift_ratio_seg_av = sum(shift_ratio,1) ./ n_segments;
    
    h_ratio_seg_av = figure;
    plot(ratio_seg_av);
    hold on
    plot(shift_ratio_seg_av,'--');
    title('Ratio');
    xlabel('Time');
    ylabel('Ratio');
    
    
    RATIO_RESULT_DIR = RESULT_DIR(1:end-1);
    while ~strcmp(RATIO_RESULT_DIR(end),filesep);
        RATIO_RESULT_DIR(end)=[];
    end
    RATIO_RESULT_DIR(end)=[];
    while ~strcmp(RATIO_RESULT_DIR(end),filesep);
        RATIO_RESULT_DIR(end)=[];
    end        

    fid_ratio = fopen([PROJECT_DIR filesep 'ratio.dat'],'a');
    if fid_ratio == -1
        error('Could not create file for ratios');
        return
    else
        fprintf(fid_ratio,'%7g  ',merg_parameters.seg_shift, shift_activity_1 shift_activity_2 shift_ratio_seg_av);
        fprintf(fid_ratio,'\n');
    end
    fclose(fid_ratio);
end












