function prPowerSpectrum(edge_parameters, merg_parameters, post_parameters, var1, var2, var1_name, var2_name)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%           Analyse Power Spectrum of time series                         %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%             Detrend data                       %%%%%%%%%%%%%%%%%%
if FILTER_METHOD > 1
    if post_parameters.do_prot
        [protrusion, protrusion_trend, av_de_trend_protrusion_time] = analyzeTimeSeries(protrusion,post_parameters.filter_method, windowSize, spline_s);
        img_protrusion  = imresize(protrusion,post_parameters.image_stretch, 'nearest');
    end
    if post_parameters.do_activity_1
        [activity_1, activity_trend_1] = analyzeTimeSeries(squeeze(activity(1,:,:)),FILTER_METHOD, windowSize, spline_s);
        [shift_activity_1, shift_activity_trend_1] = analyzeTimeSeries(squeeze(shift_activity(1,:,:)),FILTER_METHOD, windowSize, spline_s);
        %[activity, activity_trend, av_de_trend_activity_time] = analyzeTimeSeries(activity,post_parameters.filter_method, windowSize, spline_s);
        %[shift_activity, shift_activity_trend, av_de_trend_shift_activity_time] = analyzeTimeSeries(shift_activity,post_parameters.filter_method, windowSize, spline_s);
        img_activity  = imresize(activity,post_parameters.image_stretch, 'nearest');
        img_shift_protrusion  = imresize(shift_activity,post_parameters.image_stretch, 'nearest');
    end
    if post_parameters.do_retro
        [retrograde_flow, retrograde_flow_trend, av_de_trend_retrograde_flow_time] = analyzeTimeSeries(retrograde_flow,post_parameters.filter_method, windowSize, spline_s);
        [shift_retrograde_flow, shift_retrograde_flow_trend, av_de_trend_shift_retrograde_flow_time] = analyzeTimeSeries(shift_retrograde_flow,post_parameters.filter_method, windowSize, spline_s);
        img_retrograde_flow = imresize(retrograde_flow,post_parameters.image_stretch, 'nearest');
        img_shift_retrograde_flow = imresize(shift_retrograde_flow,post_parameters.image_stretch, 'nearest');
    end
    if post_parameters.do_score
        [score, score_trend, av_de_trend_score_time] = analyzeTimeSeries(score, post_parameters.filter_method, windowSize, spline_s);
        img_score = imresize(score,post_parameters.image_stretch, 'nearest');
    end

    if trend == 1
        if post_parameters.do_prot
            protrusion = protrusion_trend;
        end
        if post_parameters.do_activity_1
            activity_1 = activity_trend_1;
            shift_activity_1 = shift_activity_trend_1;
        end
        if post_parameters.do_retro
            retrograde_flow = retrograde_flow_trend;
            shift_retrograde_flow = shift_retrograde_flow_trend;
        end
        if post_parameters.do_score
            score = score_trend;
        end
    end
else
    if post_parameters.do_activity_1
        activity_1 = squeeze(activity(1,:,:));
        shift_activity_1 = squeeze(shift_activity(1,:,:));
    end
    if post_parameters.do_activity_2
        activity_2 = squeeze(activity(2,:,:));
        shift_activity_2 = squeeze(shift_activity(2,:,:));
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if 0
    % try a POD of the data!!
    [U, lambda, V] = svd(img_protrusion);

    figure,plot(diag(lambda));
    for i=1:30
        e(:,:,i) = U(:,i) * sqrt(lambda(i,i)) * (V(:,i))';
    end

    figure,imshow(e(:,:,1),[]);colormap(jet);
    figure,imshow(e(:,:,2),[]);colormap(jet);
    figure,imshow(e(:,:,3),[]);colormap(jet);

    structures  = sum(e(:,:,1:20),3);
    figure,imshow(structures,[]);colormap(jet);
    figure,imshow(img_protrusion,[]);colormap(jet);
    return;

    % n: number of realizations (time)
    % m: grid points (boundary segments)

    % 1. determine the Eigenvectors of A'A (nxn)
    % and 2. the corresponding eigenvectors V
    % protrusion = img_protrusion;
    %         AA = protrusion' * protrusion;
    %
    %         % get sigma = sqrt(lambda) of non-zero eigenvalues
    %         [U_dum, lambda_m, V] = svd(AA);
    %         lambda = diag(lambda_m,0);
    %
    %         K = U_dum(:,1) *lambda_m(1,1) * V(1,:);
    %         figure,imshow(K,[]);
    %         colormap(jet);
    %
    %         figure, imshow(protrusion,[]);
    %         colormap(jet);
    %
    %         % get the nonzero Eigenvalues
    %         [i, j, lambda_nz] = find((lambda > 10^-5) .* lambda);
    %         r = length(lambda_nz);
    %
    %         sigma = sqrt(lambda_nz);
    %         figure,plot(sigma);
    %
    %         % 3. Determine u = (1/sigma) A v
    %         for i=1:r
    %             U_under(:,i) = 1/sigma(i) .* protrusion * V(:,i);
    %         end
    %
    %         [Q,R] = qr(U_under);
    %         U_upper = Q;
    %
    %         U = [U_under U_upper];
    %
    %         r1 = U(:,1) * sigma(1) * V(1,:);
    %         r2 = U(:,2) * sigma(2) * V(2,:);
    %         r3 = U(:,3) * sigma(3) * V(3,:);
    %
    %         figure, imshow(r1);
    %         colormap(jet);
    %         figure, imshow(r2);
    %         colormap(jet);
    %         figure, imshow(r3);
    %         colormap(jet);
    %         return;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Subspace methods, also known as high-resolution methods or super-resolution methods,
% generate frequency component estimates for a signal based on an eigenanalysis or
% eigendecomposition of the correlation matrix. Examples are the multiple signal
% classification (MUSIC) method or the eigenvector (EV) method. These methods are best
% suited for line spectra - that is, spectra of sinusoidal signals - and are effective
% in the detection of sinusoids buried in noise, especially when the signal to noise
% ratios are low.

% try the music method: [S,w] = pmusic(x, p)
% each row of the input x represents a separate observation of the signal
% p is the dimension of the signal subspace. It the signal is made form one sinuosid the
% subspace is 2 because each real sinusoid is the sum of two complex
%[S,w] = pmusic(protrusion, 2);

%pmusic(x,[Inf,1.1],[],8000,7); % Window length
if post_parameters.do_prot & post_parameters.do_retro & post_parameters.do_score
    h_music_data = figure;
    subplot(3,1,1);
    pmusic(protrusion, 2)
    hold on
    subplot(3,1,2);
    pmusic(retrograde_flow, 2)
    subplot(3,1,3);
    pmusic(score, 2)

    %save data
    hgsave(h_music_data,    [RESULT_DIR 'music_data.fig']);
    print(h_music_data,     [RESULT_DIR 'music_data.eps'],'-depsc2','-tiff');
    print(h_music_data,     [RESULT_DIR 'music_data.tif'],'-dtiff');

elseif post_parameters.do_activity_1 & post_parameters.do_prot
    h_music_data = figure;
    subplot(2,1,1);
    pmusic(protrusion, 2)
    subplot(2,1,2);
    pmusic(activity_1, 2)

    %save data
    hgsave(h_music_data,    [RESULT_DIR 'music_data.fig']);
    print(h_music_data,     [RESULT_DIR 'music_data.eps'],'-depsc2','-tiff');
    print(h_music_data,     [RESULT_DIR 'music_data.tif'],'-dtiff');
elseif post_parameters.do_prot
    h_music_data = figure;
    pmusic(protrusion, 2)

    %save data
    hgsave(h_music_data,    [RESULT_DIR 'music_data.fig']);
    print(h_music_data,     [RESULT_DIR 'music_data.eps'],'-depsc2','-tiff');
    print(h_music_data,     [RESULT_DIR 'music_data.tif'],'-dtiff');
end


function F = matrix_diff(x, ref_matrix)

    [n, m] = size(ref_matrix);
    frequency=      x(1);
    average=        x(2);
    amplidude=      x(3:2+n);
    phase=          x(3+n:end);


    I=ones(n,1);
    t=linspace(0,10*(m-1),m);

    approx_matrix = average .* I * ones(1,m) + amplidude * ones(1,m).* sin(frequency*I*t + phase*ones(1,m));

    difference = ref_matrix-approx_matrix;
    F = sum(difference, 2);
