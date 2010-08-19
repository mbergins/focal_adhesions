function activity_map_interpolated = prDirInterpolateActivity(activity_map)

num_time        = size(activity_map,2);
num_segments    = size(activity_map,1);

activity_map_interpolated = activity_map;

for time = 1:num_time
    segment = 1;
    while segment <= num_segments
        if activity_map(segment,time) == -99
            segment_a = segment;
            % find next valid value
            segment = segment+1;
            while segment <= num_segments & activity_map(segment,time) == -99
                segment = segment+1;
            end
            segment_b = segment-1;
             

            % set the invalid values to the mean
            if segment_a > 1 & segment_b < num_segments
                new_score = (activity_map(segment_a-1,time) + activity_map(segment_b+1,time))/2;
                for i_seg = segment_a:segment_b
                    activity_map_interpolated(i_seg,time) = new_score;
                end
            elseif segment_a == 1 & segment_b < num_segments
                new_score = activity_map(segment_b+1,time);
                for i_seg = segment_a:segment_b
                    activity_map_interpolated(i_seg,time) = new_score;
                end
            elseif segment_a > 1 & segment_b == num_segments
                new_score = activity_map(segment_a-1,time);
                for i_seg = segment_a:segment_b
                    activity_map_interpolated(i_seg,time) = new_score;
                end
            else
                % note: you could interpolate here also
                % using the timestep before and after 
                activity_map_interpolated(:,time) = 0;
            end
        else
           segment = segment+1; 
        end % if invalid value (-99)
    end % for segments
end % for time