%source the 10th and 90th percentiles of all the adhesion intensities
min_ad_intensity = 0.2341950;
max_ad_intensity = 0.4723663;
side_ad_intensity = max_ad_intensity + 0.2;

ad_int_steps = 15;

background_mean_intensity = 0.0;

background_noise_var = 2E-3;
% background_noise_var = 0;

max_ad_size = 10;
min_ad_size = 1;

ad_padding = ceil(max_ad_size*0.4);
