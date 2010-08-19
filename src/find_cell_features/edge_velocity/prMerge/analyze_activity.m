function analyze_activity(edge_parameters, merg_parameters, post_parameters, cell_variables)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%           Analyze the data                                              %
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
DO_FILTER           = post_parameters.do_filter;
FILTER_METHOD       = post_parameters.filter_method;
trend               = post_parameters.trend;
windowSize          = post_parameters.windowSize;
spline_s            = post_parameters.spline_s;
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
end
if DO_ACTIVITY_1
   activity = cell_variables.activity;
   shift_activity = cell_variables.shift_activity;
end

x_time_axis = (merg_parameters.first_time: 1: merg_parameters.total_time_steps).* edge_parameters.time_interval;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%            Detrend   data                          %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if FILTER_METHOD > 1
    if DO_ACTIVITY_1
        [activity_1, activity_trend_1, av_de_trend_activity_time_1, av_trend_activity_time_1, av_activity_time_1] =...
                 analyzeTimeSeries(squeeze(activity(1,:,:)),FILTER_METHOD, windowSize, spline_s);
        [shift_activity_1, shift_activity_trend_1, av_de_trend_shift_activity_time_1, av_trend_shift_activity_time_1] =...
                 analyzeTimeSeries(squeeze(shift_activity(1,:,:)),FILTER_METHOD, windowSize, spline_s);

        % Plot detrended data     
        h_av_de_trend_activity_time = figure;
        plot(x_time_axis, av_de_trend_activity_time_1,'m');
        hold on
        plot(x_time_axis, av_de_trend_shift_activity_time_1,'m--');
        xlim([x_time_axis(1) x_time_axis(end)]);
        title('Detrended segment averaged activity 1');
        hgsave(h_av_de_trend_activity_time,[RESULT_DIR 'av_de_trend_activity_time_1.fig']);
        print(h_av_de_trend_activity_time, [RESULT_DIR 'av_de_trend_activity_time_1.eps'],'-depsc2','-tiff');
        print(h_av_de_trend_activity_time, [RESULT_DIR 'av_de_trend_activity_time_1.tif'],'-dtiff');
        
        % Plot trend
        h_av_trend_activity_time = figure;
        plot(x_time_axis, av_trend_activity_time_1,'m');
        hold on
        %plot(x_time_axis, av_activity_time_1./av_trend_activity_time_1(1),'m');
        plot(x_time_axis, av_trend_shift_activity_time_1,'m--');
        xlim([x_time_axis(1) x_time_axis(end)]);
        title('Trend segment averaged activity 1'); 
        hgsave(h_av_trend_activity_time,[RESULT_DIR 'av_trend_activity_time_1.fig']);
        print(h_av_trend_activity_time, [RESULT_DIR 'av_trend_activity_time_1.eps'],'-depsc2','-tiff');
        print(h_av_trend_activity_time, [RESULT_DIR 'av_trend_activity_time_1.tif'],'-dtiff');
    end

    if DO_SCORE
        [score, score_trend, av_de_trend_score_time, av_trend_score_time, av_score_time, var_score_trend, var_score_detrend] =...
                                                        analyzeTimeSeries(score,FILTER_METHOD, windowSize, spline_s);
                                                    
        [shift_score, shift_score_trend, av_de_trend_shift_score_time, av_trend_shift_score_time, av_shift_score_time, var_shift_score_trend, var_shift_score_detrend] =...
                                                        analyzeTimeSeries(shift_score,FILTER_METHOD, windowSize, spline_s);                                                        

        av_de_trend_score_time = nanmean(score,1);
        h_av_de_trend_score_time = figure;
        plot(x_time_axis, av_de_trend_score_time,'b');
        xlim([x_time_axis(1) x_time_axis(end)]);
        title('De-trended segment averaged scores');
        hgsave(h_av_de_trend_score_time,[RESULT_DIR 'av_de_trend_score_time.fig']);
        print(h_av_de_trend_score_time, [RESULT_DIR 'av_de_trend_score_time.eps'],'-depsc2','-tiff');
        print(h_av_de_trend_score_time, [RESULT_DIR 'av_de_trend_score_time.tif'],'-dtiff');
        
        % plot trend
        h_av_trend_score_time = figure;
        plot(x_time_axis, av_trend_score_time,'b--');
        hold on
        plot(x_time_axis, av_score_time,'b');
        xlim([x_time_axis(1) x_time_axis(end)]);
        title('Trend segment averaged scores');
        legend('Trend', 'Original');
        hgsave(h_av_trend_score_time,[RESULT_DIR 'av_trend_score_time.fig']);
        print(h_av_trend_score_time, [RESULT_DIR 'av_trend_score_time.eps'],'-depsc2','-tiff');
        print(h_av_trend_score_time, [RESULT_DIR 'av_trend_score_time.tif'],'-dtiff');         
    end

    if DO_RETRO
        [retrograde_flow, retrograde_flow_trend, av_de_trend_retrograde_flow_time,...
         av_trend_retrograde_flow_time, av_retrograde_flow_time,...
         var_retrograde_flow_trend, var_retrograde_flow_detrend] = ...
               analyzeTimeSeries(retrograde_flow,FILTER_METHOD, windowSize, spline_s);
           
        [shift_retrograde_flow, shift_retrograde_flow_trend, av_de_trend_shift_retrograde_flow_time] = analyzeTimeSeries(shift_retrograde_flow,FILTER_METHOD, windowSize, spline_s);

        av_de_trend_retrograde_flow_time = sum(retrograde_flow,1) ./ sum(retrograde_flow ~=0, 1);
        av_de_trend_shift_retrograde_flow_time = sum(shift_retrograde_flow,1) ./ sum(shift_retrograde_flow ~=0, 1);

        h_av_de_trend_retro_time = figure;
        plot(x_time_axis, av_de_trend_retrograde_flow_time,'r');
        hold on
        plot(x_time_axis, av_de_trend_shift_retrograde_flow_time,'r--');
        xlim([x_time_axis(1) x_time_axis(end)]);
        legend('retro','shift flow');
        title('Detrended segment averaged retrograde flow');
        hgsave(h_av_de_trend_retro_time,[RESULT_DIR 'av_de_trend_retro_time.fig']);
        print(h_av_de_trend_retro_time, [RESULT_DIR 'av_de_trend_retro_time.eps'],'-depsc2','-tiff');
        print(h_av_de_trend_retro_time, [RESULT_DIR 'av_de_trend_retro_time.tif'],'-dtiff');
        
        % plot trend
        h_av_trend_retro_time = figure;
        plot(x_time_axis, av_trend_retrograde_flow_time,'r--');
        hold on
        plot(x_time_axis, av_retrograde_flow_time,'r');
        xlim([x_time_axis(1) x_time_axis(end)]);
        title('Trend segment averaged retrograde flow');
        hgsave(h_av_trend_retro_time,[RESULT_DIR 'av_trend_retro_time.fig']);
        print(h_av_trend_retro_time, [RESULT_DIR 'av_trend_retro_time.eps'],'-depsc2','-tiff');
        print(h_av_trend_retro_time, [RESULT_DIR 'av_trend_retro_time.tif'],'-dtiff');     
    end

    if DO_PROT
        [protrusion, protrusion_trend, av_de_trend_protrusion_time,...
            av_trend_protrusion_time, av_protrusion_time, var_protrusion_trend, var_protrusion_detrend] =...
            analyzeTimeSeries(protrusion,FILTER_METHOD, windowSize, spline_s);

        % av_de_trend_protrusion_time = sum(protrusion,1) ./ sum(protrusion ~=0, 1);
        h_av_de_trend_protrusion_time = figure;
        plot(x_time_axis, av_de_trend_protrusion_time,'g');
        xlim([x_time_axis(1) x_time_axis(end)]);
        title('Detrended segment averaged protrusion');
        hgsave(h_av_de_trend_protrusion_time,[RESULT_DIR 'av_de_trend_protrusion_time.fig']);
        print(h_av_de_trend_protrusion_time, [RESULT_DIR 'av_de_trend_protrusion_time.eps'],'-depsc2','-tiff');
        print(h_av_de_trend_protrusion_time, [RESULT_DIR 'av_de_trend_protrusion_time.tif'],'-dtiff');
        
        % plot trend
        h_av_trend_protrusion_time = figure;
        plot(x_time_axis, av_trend_protrusion_time,'g--');
        hold on
        plot(x_time_axis, av_protrusion_time,'g');
        xlim([x_time_axis(1) x_time_axis(end)]);
        title('Trend segment averaged protrusion');
        hgsave(h_av_trend_protrusion_time,[RESULT_DIR 'av_trend_protrusion_time.fig']);
        print(h_av_trend_protrusion_time, [RESULT_DIR 'av_trend_protrusion_time.eps'],'-depsc2','-tiff');
        print(h_av_trend_protrusion_time, [RESULT_DIR 'av_trend_protrusion_time.tif'],'-dtiff');        
    end
    
    if DO_ACTIVITY_1 & DO_PROT
        % Compare trend
        h_av_trend_time = figure;
        plot(x_time_axis, av_trend_protrusion_time,'g');
        xlim([x_time_axis(1) x_time_axis(end)]);
        hold on
        %plot(x_time_axis, av_protrusion_time,'g--');
        % Set the second y axis for activity
        ax1 = gca;
        set(ax1,'box','off');
        xlabel('Time [s]');
        ylabel('Protrusion [nm/s]');
        ax2 = axes('position',get(ax1,'position'));
        set(ax2,'XAxisLocation','bottom','YAxisLocation','right','color','no');
        
        % plot activity
        line(x_time_axis, av_trend_activity_time_1,'Color','m');
        %line(x_time_axis, av_trend_shift_activity_time_1,'Color','m','LineStyle','--');
        xlim([x_time_axis(1) x_time_axis(end)]);
        %line(x_time_axis, av_activity_time,'Color','m','LineStyle','-');
        ylabel('Activity 1 [AU]');
        hgsave(h_av_trend_time,[RESULT_DIR 'av_trend_time.fig']);
        print(h_av_trend_time, [RESULT_DIR 'av_trend_time.eps'],'-depsc2','-tiff');
        print(h_av_trend_time, [RESULT_DIR 'av_trend_time.tif'],'-dtiff');      
        
        % Compare detrended data
        h_av_de_trend_prot_activity_time = figure;
        plot(x_time_axis, av_de_trend_activity_time_1,'m');
        hold on
        %plot(x_time_axis, av_de_trend_shift_activity_time_1,'--m');
        plot(x_time_axis, av_de_trend_protrusion_time,'g');
        xlim([x_time_axis(1) x_time_axis(end)]);
        title('Detrended segment averaged protrusion and activity 1');
        hgsave(h_av_de_trend_prot_activity_time,[RESULT_DIR 'av_de_trend_prot_activity_time_1.fig']);
        print(h_av_de_trend_prot_activity_time, [RESULT_DIR 'av_de_trend_prot_activity_time_1.eps'],'-depsc2','-tiff');
        print(h_av_de_trend_prot_activity_time, [RESULT_DIR 'av_de_trend_aprot_ctivity_time_1.tif'],'-dtiff');  
    end
    
    if DO_PROT & DO_SCORE & DO_RETRO 
        % plot normalized trend
        %prestd 
        % Preprocess data so that its mean is 0 and the standard deviation is 1
        h_av_trend_comp = figure;
        plot(x_time_axis, prestd(av_trend_score_time),'b');
        hold on
        plot(x_time_axis, prestd(av_trend_retrograde_flow_time),'r');
        plot(x_time_axis, prestd(av_trend_protrusion_time),'g');
        title('Comparison of the trend');
        hgsave(h_av_trend_comp,[RESULT_DIR 'av_trend_comp.fig']);
        print(h_av_trend_comp, [RESULT_DIR 'av_trend_comp.eps'],'-depsc2','-tiff');
        print(h_av_trend_comp, [RESULT_DIR 'av_trend_comp.tif'],'-dtiff');       
        
        
        h_av_detrend_comp = figure;
        plot(x_time_axis, prestd(av_de_trend_score_time),'b');
        hold on
        plot(x_time_axis, prestd(av_de_trend_retrograde_flow_time),'r');
        plot(x_time_axis, prestd(av_de_trend_protrusion_time),'g');
        title('Comparison of the detrend');
        hgsave(h_av_detrend_comp,[RESULT_DIR 'av_detrend_comp.fig']);
        print(h_av_detrend_comp, [RESULT_DIR 'av_detrend_comp.eps'],'-depsc2','-tiff');
        print(h_av_detrend_comp, [RESULT_DIR 'av_detrend_comp.tif'],'-dtiff');       
        
        % check the alpha factor
        alpha_trend = prestd(av_trend_protrusion_time) - prestd(av_trend_retrograde_flow_time);
        h_alpha_trend = figure;
        plot(x_time_axis, alpha_trend,'c');
        hold on
        plot(x_time_axis, prestd(av_de_trend_score_time),'b');
        title('Alpha factor');
    end
    
    
    if DO_PROT & DO_SCORE
        % plot normalized trend
        %prestd 
        % Preprocess data so that its mean is 0 and the standard deviation is 1
        h_av_trend_comp2 = figure;
        plot(x_time_axis, prestd(av_trend_protrusion_time),'g');
        hold on
        plot(x_time_axis, prestd(av_trend_score_time),'b');
        plot(x_time_axis, prestd(av_trend_shift_score_time),'--b');
        title('Comparison of the trend');
        legend('Protrusion', 'Scores', 'Shifted scores');
        hgsave(h_av_trend_comp2,[RESULT_DIR 'av_trend_comp.fig']);
   
        h_av_detrend_comp2 = figure;
        plot(x_time_axis, prestd(av_de_trend_score_time),'b');
        hold on
        plot(x_time_axis, prestd(av_de_trend_protrusion_time),'g');
        title('Comparison of the de-trend');
        hgsave(h_av_detrend_comp2,[RESULT_DIR 'av_detrend_comp.fig']);      
    end    
    
    if DO_PROT & DO_RETRO
        % plot normalized trend
        %prestd 
        % Preprocess data so that its mean is 0 and the standard deviation is 1
        h_av_trend_comp3 = figure;
        plot(x_time_axis, prestd(av_trend_protrusion_time),'g');
        hold on
        plot(x_time_axis, prestd(av_trend_retrograde_flow_time),'b');
        %plot(x_time_axis, prestd(av_trend_shift_score_time),'--b');
        title('Comparison of the trend');
        legend('Protrusion', 'Flow');
        hgsave(h_av_trend_comp3,[RESULT_DIR 'av_trend_comp3.fig']);
   
%         h_av_detrend_comp3 = figure;
%         plot(x_time_axis, prestd(av_de_trend_score_time),'b');
%         hold on
%         plot(x_time_axis, prestd(av_de_trend_protrusion_time),'g');
%         title('Comparison of the de-trend');
%         hgsave(h_av_detrend_comp2,[RESULT_DIR 'av_detrend_comp3.fig']);      
    end      
    
end % detrend

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
