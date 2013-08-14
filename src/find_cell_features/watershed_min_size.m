function final_label_mat = watershed_min_size(image,label_mat,min_size)

%This function takes advantage of the fact that with a minimum size
%threshold, the watershed segmentation only needs to consider objects at
%least 2 times the minimum size threshold for the segmentation
props = regionprops(label_mat,image,'Area');
area_vals = [props.Area];
large_ad_nums = find(area_vals >= 2*min_size);

final_label_mat = zeros(size(label_mat,1),size(label_mat,2));
for this_ad_num = 1:max(label_mat(:))
    this_ad = label_mat == this_ad_num;
    
    if (mod(this_ad_num,round(max(label_mat(:))/10)) == 0) 
        fprintf('Done with %d/%d adhesions\n',this_ad_num,max(label_mat(:)));
    end
    
    %we aren't dealing with a large adhesion, so it can't be split by the
    %watershed methods, add it to the final label matrix and move on
    if (not(any(this_ad_num == large_ad_nums)))
        final_label_mat(this_ad) = max(final_label_mat(:)) + 1;
        continue;
    end
    
    %we need to determine the order in which to add the pixels, brightest
    %pixels will be in first, let's find the linear indexes and put them in
    %order
    sorted_ad_intensities = sort(image(this_ad),'descend');
    
    %this for loop needs to deal with the unlikely case of ties in ad
    %intensities, in that case, take the first one out of the find and set
    %it to an impossible value, so it won't match on the next intensity
    %search
    sorted_linear_indexes = zeros(size(sorted_ad_intensities));
    this_ad_image = this_ad .* image;
    for ad_int_ind = 1:length(sorted_ad_intensities)
        sorted_linear_indexes(ad_int_ind) = find(this_ad_image == sorted_ad_intensities(ad_int_ind),1,'first');        
        this_ad_image(sorted_linear_indexes(ad_int_ind)) = -Inf;
    end
    
    %Up till the minimum independent size barrier, we can't have any
    %watershed splitting events, so go ahead and seed in up to the minimum
    %size pixels and label them using bwlabel
    watershed_ad = zeros(size(label_mat,1),size(label_mat,2));
    for i=1:min_size
       watershed_ad(sorted_linear_indexes(1)) = 1;
       sorted_linear_indexes = sorted_linear_indexes(2:end);
    end
    
    watershed_ad = bwlabel(watershed_ad,4);
    
    while (length(sorted_linear_indexes) >= 1)
        watershed_ad = add_single_pixel_watershed(watershed_ad,sorted_linear_indexes(1),min_size);
        sorted_linear_indexes = sorted_linear_indexes(2:end);
    end
    
    %add the watershed segmented FA back to the master segmentation
    unique_labels = unique(watershed_ad(:));
    for j=2:length(unique_labels);
        final_label_mat(watershed_ad == unique_labels(j)) = max(final_label_mat(:)) + 1;
    end
end


function watershed_ads = add_single_pixel_watershed(watershed_ads,pix_pos,min_size)

initial_size = sum(sum(watershed_ads > 0));

[pix_pos_ind(1),pix_pos_ind(2)] = ind2sub(size(watershed_ads),pix_pos);

touching_ad_nums = zeros(4,1);
%wrap these calls in try, since they could attempt to access indexes
%outside the 1-size range
try touching_ad_nums(1) = watershed_ads(pix_pos_ind(1) - 1,pix_pos_ind(2)); end %#ok<TRYNC>
try touching_ad_nums(2) = watershed_ads(pix_pos_ind(1) + 1,pix_pos_ind(2)); end %#ok<TRYNC>
try touching_ad_nums(3) = watershed_ads(pix_pos_ind(1),pix_pos_ind(2) - 1); end %#ok<TRYNC>
try touching_ad_nums(4) = watershed_ads(pix_pos_ind(1),pix_pos_ind(2) + 1); end %#ok<TRYNC>

touching_ad_nums = nonzeros(unique(touching_ad_nums));
%if there aren't any pixels which were connected to newest pixel, add the
%newest pixel as a new adhesion, otherwise, start a more complicated
%procedure
if (sum(touching_ad_nums) == 0)
    watershed_ads(pix_pos) = max(watershed_ads(:)) + 1;
else
    %if there is only one touching adhesion, label the new pixel with that
    %adhesion's label
    if (length(touching_ad_nums) == 1)
        watershed_ads(pix_pos) = touching_ad_nums(1);
    else
        touching_ads = ismember(watershed_ads,touching_ad_nums).*watershed_ads;
        
        props = regionprops(touching_ads,'Area','Centroid');
        
        %Check if all the connected adhesions are below the minimum size,
        %if they are, merge all the adhesions, otherwise, trigger a more
        %complicated procedure
        meets_size_min = [props.Area] >= min_size;
        
        if (all([props.Area] < min_size))
            min_ad_num = min(touching_ad_nums);
            watershed_ads(touching_ads > 0) = min_ad_num;
            watershed_ads(pix_pos) = min_ad_num;
            1;
        elseif (sum(meets_size_min) == 1)
            %There is only one large adhesion touching, find the adhesion
            %number and label the new pixel with that number
            large_ad_num = find(meets_size_min);
            watershed_ads(pix_pos) = large_ad_num;
            
            doesnt_meet_size_min = find([props.Area] > 0 & [props.Area] < min_size);
            if (not(isempty(doesnt_meet_size_min)))
                %now assign those touching positions to the same adhesion
                %as selected by the distance metric
                for i=doesnt_meet_size_min
                    watershed_ads(watershed_ads == i) = large_ad_num;
                end
            end
        else
            %we now have at least two large adhesions touching the current
            %pixel, we will decide on the assignment based on the distance
            %to the centroids of the current objects
            centroid_pos = reshape([props.Centroid],2,[]);
            dists = (centroid_pos(1,:) - pix_pos_ind(2)).^2 + (centroid_pos(2,:) - pix_pos_ind(1)).^2;
            
            %we need to block any small adhesion distances that might also
            %be touching the current position, highly unlikely, but just to
            %make sure...
            dists = dists.*meets_size_min;
            
            min_dist_ad_num = find(dists == min(dists),1,'first');
            watershed_ads(pix_pos) = min_dist_ad_num;
            
            %now do a search for any small adhesions that didn't meet the
            %minimum size threshold, but are still touching the current
            %position, again, highly unlikely, but just in case...
            doesnt_meet_size_min = find([props.Area] > 0 & [props.Area] < min_size);
            
            if (not(isempty(doesnt_meet_size_min)))
                %now assign those touching positions to the same adhesion
                %as selected by the distance metric
                for i=doesnt_meet_size_min
                    watershed_ads(watershed_ads == i) = min_dist_ad_num;
                end
            end
        end

    end
end

final_size = sum(sum(watershed_ads > 0));
assert(initial_size + 1 == final_size)
