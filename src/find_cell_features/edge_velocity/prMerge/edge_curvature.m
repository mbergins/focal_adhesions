function edge_curvature(edge_parameters, merg_parameters, post_parameters, cell_variables)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%           Analyze data using FFT algorithm                              %
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

delta_t = edge_parameters.time_interval;




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% Get the curvature of the splines in each time step %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if DO_ACTIVITY_1
    h_fft_activity = figure;
    subplot(2,1,1)
    title('Activity 1');
    plot(frequ(2:floor(end/2)), av_activity_fft_mag_1(2:floor(end/2)),'m');
    subplot(2,1,2)
    plot(frequ(2:floor(end/2)), av_activity_fft_phase_1(2:floor(end/2)),'m');
    xlabel('Frequency [1/sec]');
    ylabel('[-]');
    hgsave(h_fft_activity,   [RESULT_DIR 'av_fft_activity_1.fig']);
    print(h_fft_activity,    [RESULT_DIR 'av_fft_activity_1.eps'],'-depsc2','-tiff');
    print(h_fft_activity,    [RESULT_DIR 'av_fft_activity_1.tif'],'-dtiff');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%













%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%   Matlab spectrum estimator %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if 0
    protrusion_periodogram              = zeros(fft_length+1,1);
    retrograde_flow_periodogram         = zeros(fft_length+1,1);
    shift_retrograde_flow_periodogram   = zeros(fft_length+1,1);
    scores_periodogram                  = zeros(fft_length+1,1);

    for i=1:size(protrusion,1)
        %[temp_perio, w_w]  = pwelch(protrusion(i,:));
        %         [temp1, w_p]    = periodogram(protrusion(i,:));
        %         [temp2, w_p]    = periodogram(retrograde_flow(i,:));
        %         [temp3, w_p]    = periodogram(shift_retrograde_flow(i,:));
        %         [temp4, w_p]    = periodogram(scores(i,:));
        %         [temp1, w_p]    = pburg(protrusion(i,:),30);
        %         [temp2, w_p]    = pburg(retrograde_flow(i,:),30);
        %         [temp3, w_p]    = pburg(shift_retrograde_flow(i,:),30);
        %         [temp4, w_p]    = pburg(scores(i,:),30);

        mod_order = 6;

        [temp1, w_p]    = pmusic(protrusion(i,:),mod_order, 2*fft_length, 1/delta_t);
        [temp2, w_p]    = pmusic(retrograde_flow(i,:),mod_order, 2*fft_length, 1/delta_t);
        [temp3, w_p]    = pmusic(shift_retrograde_flow(i,:),mod_order, 2*fft_length, 1/delta_t);
        [temp4, w_p]    = pmusic(scores(i,:),mod_order, 2*fft_length, 1/delta_t);

        protrusion_periodogram              = protrusion_periodogram            + temp1;
        retrograde_flow_periodogram         = retrograde_flow_periodogram       + temp2;
        shift_retrograde_flow_periodogram   = shift_retrograde_flow_periodogram + temp3;
        scores_periodogram                  = scores_periodogram                + temp4;
    end

    protrusion_periodogram              = protrusion_periodogram            ./ size(protrusion,1);
    retrograde_flow_periodogram         = retrograde_flow_periodogram       ./ size(protrusion,1);
    shift_retrograde_flow_periodogram   = shift_retrograde_flow_periodogram ./ size(protrusion,1);
    scores_periodogram                  = scores_periodogram                ./ size(protrusion,1);

    h_spec_est = figure;
    plot(w_p(2:end),protrusion_periodogram(2:end),'g');
    hold on
    plot(w_p(2:end),retrograde_flow_periodogram(2:end),'r');
    plot(w_p(2:end),shift_retrograde_flow_periodogram(2:end),'c');
    plot(w_p(2:end),scores_periodogram(2:end),'b');
    title('Spectrum estimation');
    xlabel('Frequency [1/sec]');
    ylabel('[-]');
    legend('Protrusion', 'Retrograde flow', 'Shift retrograde flow', 'Scores');

    hgsave(h_spec_est,   [RESULT_DIR 'spec_est.fig']);
    print(h_spec_est,    [RESULT_DIR 'spec_est.eps'],'-depsc2','-tiff');
    print(h_spec_est,    [RESULT_DIR 'spec_est.tif'],'-dtiff');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% here we remodel the matrix into a array and try to this way a FFT
%     protrusion_array = reshape(protrusion,[1 prod(size(protrusion))]);
%
%     protrusion_fft                  = fft(protrusion_array);
%     av_protrusion_fft_mag           = abs(protrusion_fft);
%     av_protrusion_fft_phase         = unwrap(angle(protrusion_fft));
%     frequ = (1:floor(size(protrusion_fft,2)/2)) ./
%     (size(protrusion_array,2) * prot_parameter.time_interval);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%