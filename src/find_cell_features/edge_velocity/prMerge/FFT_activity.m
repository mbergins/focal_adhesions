function FFT_activity(edge_parameters, merg_parameters, post_parameters, cell_variables)
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
%%%%%             Detrend data                       %%%%%%%%%%%%%%%%%%
if TREND_METHOD > 1
    if DO_PROT
        [protrusion, protrusion_trend, av_de_trend_protrusion_time] = analyzeTimeSeries(protrusion,TREND_METHOD, windowSize, spline_s);
    end
    if DO_ACTIVITY_1
        [activity_1, activity_trend_1, av_de_trend_activity_time_1] = analyzeTimeSeries(squeeze(activity(1,:,:)),TREND_METHOD, windowSize, spline_s);
        [shift_activity_1, shift_activity_trend_1, av_de_trend_shift_activity_time_1] = analyzeTimeSeries(squeeze(shift_activity(1,:,:)),TREND_METHOD, windowSize, spline_s);
    end
    if DO_RETRO
        [retrograde_flow, retrograde_flow_trend, av_de_trend_retrograde_flow_time] = analyzeTimeSeries(retrograde_flow,TREND_METHOD, windowSize, spline_s);
        [shift_retrograde_flow, shift_retrograde_flow_trend, av_de_trend_shift_retrograde_flow_time] = analyzeTimeSeries(shift_retrograde_flow,TREND_METHOD, windowSize, spline_s);
    end
    if DO_SCORE
        [score, score_trend, av_de_trend_score_time] = analyzeTimeSeries(score, TREND_METHOD, windowSize, spline_s);
        [shift_score, shift_score_trend, av_de_trend_shift_score_time] = analyzeTimeSeries(shift_score, TREND_METHOD, windowSize, spline_s);
    end
    
    if trend == 1
        if DO_PROT
            protrusion = protrusion_trend;
        end
        if DO_ACTIVITY_1
            activity_1 = activity_trend_1;
            shift_activity_1 = shift_activity_trend_1;
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
    end        
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% Do the FFT analysis directly on the data   %%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate the mean values
% equally you could remove the first component of the FFT
% if there is a trend in the data, you should use a regression
% or a polynomial approximation
% mean_prot = abs(robust_av_protrusion(1)+0.5*robust_av_protrusion(2)*x_time_axis(end)) + 3*robust_sigma_protrusion;
%mean_prot                   = mean(protrusion,2);
%mean_retrograde_flow        = mean(retrograde_flow,2);
%mean_shift_retrograde_flow  = mean(shift_retrograde_flow,2);
%mean_scores                 = mean(scores,2);


% Find the next power to our array langth
% 64, 128, 512,1024,2048
fft_length = 2^nextpow2(protrusion(1,:));

% Construct a window
%w1 = hann(size(protrusion,2))';
w1 = kaiser(size(protrusion,2),4)';
w1 = ones(1,size(protrusion,2));
w  = repmat(w1,size(protrusion,1),1);

%     x = 1:size(protrusion,2);
%     x=x.*5;
%     y = cos(x*2*pi/30) + cos(x*2*pi/60+27);
%     figure
%     plot(x,y);
%     protrusion=repmat(y,size(protrusion,1),1);

if DO_PROT
    protrusion_fft                  = fft(w .* protrusion,fft_length,2);
    av_protrusion_fft               = mean(protrusion_fft, 1);
    av_protrusion_fft_mag           = abs(av_protrusion_fft);
    av_protrusion_fft_phase         = unwrap(angle(av_protrusion_fft));
end
if DO_RETRO
    retrograde_flow_fft             = fft(w .* retrograde_flow,fft_length,2);
    av_retrograde_flow_fft          = mean(retrograde_flow_fft, 1);
    av_retrograde_flow_fft_mag      = abs(av_retrograde_flow_fft);
    av_retrograde_flow_fft_phase    = unwrap(angle(av_retrograde_flow_fft));%

    shift_retrograde_flow_fft       = fft(w .* shift_retrograde_flow,fft_length,2);
    av_shift_retrograde_flow_fft    = mean(shift_retrograde_flow_fft, 1);
    av_shift_retrograde_flow_fft_mag= abs(av_shift_retrograde_flow_fft);
    av_shift_retrograde_flow_fft_phase= unwrap(angle(av_shift_retrograde_flow_fft));%
end
if DO_SCORE
    scores_fft                      = fft(w .* score,fft_length,2);
    av_scores_fft                   = mean(scores_fft, 1);
    av_scores_fft_mag               = abs(av_scores_fft);
    av_scores_fft_phase             = unwrap(angle(av_scores_fft));
    
    shift_scores_fft                      = fft(w .* shift_score,fft_length,2);
    av_shift_scores_fft                   = mean(shift_scores_fft, 1);
    av_shift_scores_fft_mag               = abs(av_shift_scores_fft);
    av_shift_scores_fft_phase             = unwrap(angle(av_shift_scores_fft));    
end
if DO_ACTIVITY_1
    activity_fft_1                      = fft(w .* activity_1,fft_length,2);
    av_activity_fft_1                   = mean(activity_fft_1, 1);
    av_activity_fft_mag_1               = abs(av_activity_fft_1);
    av_activity_fft_phase_1             = unwrap(angle(av_activity_fft_1));
end

frequ = (1:fft_length) ./ (fft_length * delta_t);


if DO_PROT
    h_fft_prot = figure;
    subplot(2,1,1)
    title('Protrusion');
    plot(frequ(2:floor(end/2)), av_protrusion_fft_mag(2:floor(end/2)),'g');
    subplot(2,1,2)
    plot(frequ(2:floor(end/2)), av_protrusion_fft_phase(2:floor(end/2)),'g');
    hgsave(h_fft_prot,   [RESULT_DIR 'av_fft_prot.fig']);
    print(h_fft_prot,    [RESULT_DIR 'av_fft_prot.eps'],'-depsc2','-tiff');
    print(h_fft_prot,    [RESULT_DIR 'av_fft_prot.tif'],'-dtiff');
end

if DO_RETRO
    h_fft_retro = figure;
    subplot(2,1,1)
    title('Retrograde flow');
    plot(frequ(2:floor(end/2)), av_retrograde_flow_fft_mag(2:floor(end/2)),'r');
    plot(frequ(2:floor(end/2)), av_shift_retrograde_flow_fft_mag(2:floor(end/2)),'r');
    subplot(2,1,2)
    plot(frequ(2:floor(end/2)), av_retrograde_flow_fft_phase(2:floor(end/2)),'r');
    plot(frequ(2:floor(end/2)), av_shift_retrograde_flow_fft_phase(2:floor(end/2)),'r');
    hgsave(h_fft_retro,   [RESULT_DIR 'av_fft_retro.fig']);
    print(h_fft_retro,    [RESULT_DIR 'av_fft_retro.eps'],'-depsc2','-tiff');
    print(h_fft_retro,    [RESULT_DIR 'av_fft_retro.tif'],'-dtiff');
end

if DO_SCORE
    h_fft_score = figure;
    subplot(2,1,1)
    title('Polymerization');
    plot(frequ(2:floor(end/2)), av_scores_fft_mag(2:floor(end/2)),'b');
    subplot(2,1,2)
    plot(frequ(2:floor(end/2)), av_scores_fft_phase(2:floor(end/2)),'b');
    hgsave(h_fft_score,   [RESULT_DIR 'av_fft_score.fig']);
    print(h_fft_score,    [RESULT_DIR 'av_fft_score.eps'],'-depsc2','-tiff');
    print(h_fft_score,    [RESULT_DIR 'av_fft_score.tif'],'-dtiff');
end

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