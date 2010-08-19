function coupling_activity(edge_parameters, merg_parameters, post_parameters, cell_variables)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%           Estimate the auto correlations and the                        %
%           cross correlations                                            %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PR_ALPHA_DIR        = post_parameters.pr_alpha_dir;  
PROJECT_DIR         = post_parameters.project_dir;  
POST_PROC_OPTIONS   = post_parameters.post_proc_options;  
DO_PROT             = post_parameters.do_prot;  
DO_RETRO            = post_parameters.do_retro;  
DO_SCORE            = post_parameters.do_score;  
DO_ACTIVITY_1       = post_parameters.do_activity_1;  
DO_ACTIVITY_2       = post_parameters.do_activity_2;  
DO_DETREND          = post_parameters.do_detrend;
TREND_METHOD        = post_parameters.trend_method;
trend               = post_parameters.trend;
windowSize          = post_parameters.windowSize;
spline_s            = post_parameters.spline_s;
corr_lag            = post_parameters.corr_lag;
RESULT_DIR          = post_parameters.result_dir;


if DO_PROT
    protrusion = cell_variables.protrusion;
end
if DO_SCORE
    score = cell_variables.score;
    shift_score = cell_variables.shift_score;
end
if DO_RETRO
    retrograde_flow = cell_variables.retrograde_flow;
    shift_retrograde_flow = cell_variables.shift_retrograde_flow;
    if isfield(cell_variables,'vector_field')
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%  Lin here is you data  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % the format is (seg_index,time_index,coordinate)
        vector_field = cell_variables.vector_field; 
        shift_vector_field = cell_variables.shift_vector_field;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    else
        disp('No vector data');
        return;
    end
end
if DO_ACTIVITY_1 | DO_ACTIVITY_2
   activity = cell_variables.activity;
   shift_activity = cell_variables.shift_activity;
end

% get the number of segments
% get the number of time steps
if DO_PROT
    n_segments = size(protrusion, 1);
    n_time     = size(protrusion, 2);
else
    n_segments = size(activity, 2);
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
        [shift_score, shift_score_trend, av_de_trend_shift_score_time] = analyzeTimeSeries(shift_score,TREND_METHOD, windowSize, spline_s);
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
            shift_score = shift_score_trend;
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

% Calculate the significance bounds
% 99% probability -> 2.58
% 95% probability -> 1.96
% 90% probability -> 1.645

upper_sig_bound = (n_time)^(-0.5) * 1.96;
lower_sig_bound = -upper_sig_bound;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%  Put your code here %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%