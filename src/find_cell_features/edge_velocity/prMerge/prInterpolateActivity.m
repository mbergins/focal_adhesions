function activity_map_interpolated = prInterpolateActivity(activity_map)

i=1;
num_time        = size(activity_map,2);
num_segments    = size(activity_map,1);

for time = 1:num_time
    for segment = 1:num_segments
        if activity_map(segment,time) ~= -99
            activity_map_vec(i) = activity_map(segment,time);
            x(i) = segment;
            y(i) = time;
            i=i+1;
        end
    end
end
[XI,YI] = meshgrid(1:num_segments, 1:num_time);
activity_map_interpolated = griddata(x,y,activity_map_vec ,XI,YI,'v4');
activity_map_interpolated = activity_map_interpolated';

% filter out NaN
[x_nan,y_nan]=find(isnan(activity_map));
for i=1:length(x_nan)
    activity_map_interpolated(x_nan(i), y_nan(i)) = 0;
end

% figure
% surface(XI,YI,activity_map_interpolated');
% hold on
% plot3(x,y,activity_map_vec,'o');