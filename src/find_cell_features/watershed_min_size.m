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
    
    watershed_ad = zeros(size(label_mat,1),size(label_mat,2));
    while (length(sorted_linear_indexes) >= 1)
        watershed_ad = add_single_pixel(watershed_ad,sorted_linear_indexes(1),min_size);
        sorted_linear_indexes = sorted_linear_indexes(2:end);
    end
    
    %add the watershed segmented FA back to the master segmentation
    unique_labels = unique(watershed_ad(:));
    assert(unique_labels(1) == 0)
    for j=2:length(unique_labels);
        final_label_mat(watershed_ad == unique_labels(j)) = max(final_label_mat(:)) + 1;
    end
end

function ad_zamir = add_single_pixel(ad_zamir,pix_pos,min_size)

[pix_pos_ind(1),pix_pos_ind(2)] = ind2sub(size(ad_zamir),pix_pos);

%save the number of pixels in currently in adhesions, this will be compared
%at the end of processing to make sure another pixel has been added
initial_size = sum(sum(im2bw(ad_zamir,0)));

%now locate the adhesions in the current adhesions that touch the newest
%selected pixel, store these adhesions in connected_ad
connected_ad = false(size(ad_zamir));
connected_ad(pix_pos) = 1;
ad_nums = zeros(4);
%wrap these calls in try, since they could attempt to access indexes
%outside the 1-size range
try ad_nums(1) = ad_zamir(pix_pos_ind(1) - 1,pix_pos_ind(2)); end %#ok<TRYNC>
try ad_nums(2) = ad_zamir(pix_pos_ind(1) + 1,pix_pos_ind(2)); end %#ok<TRYNC>
try ad_nums(3) = ad_zamir(pix_pos_ind(1),pix_pos_ind(2) - 1); end %#ok<TRYNC>
try ad_nums(4) = ad_zamir(pix_pos_ind(1),pix_pos_ind(2) + 1); end %#ok<TRYNC>

for i = 1:length(ad_nums)
    if (ad_nums(i) ~= 0)
        connected_ad(ad_zamir == ad_nums(i)) = 1;
    end
end

%build a binary image of the current touching adhesions, without the new
%pixel
old_ad = false(size(ad_zamir));
old_ad(and(ad_zamir > 0,connected_ad)) = 1;

%if there aren't any pixels which were connected to newest pixel, add the
%newest pixel as a new adhesion, otherwise, start a more complicated
%procedure
if (sum(old_ad(:)) == 0)
    ad_zamir(connected_ad) = max(ad_zamir(:)) + 1;
else
    %build an image with the connected old adhesions renumbered to start
    %with one, since it will be easier to work with the regionprops output when
    %numbered from one
    relabeled_old_ad = ad_zamir;
    relabeled_old_ad(old_ad ~= 1) = 0;
    ad_nums = unique(relabeled_old_ad);
    assert(ad_nums(1) == 0, 'Error in collecting relabeled_old_ad unique ad numbers')
    for i = 2:length(ad_nums)
        relabeled_old_ad(relabeled_old_ad == ad_nums(i)) = i - 1;
    end
    assert(all(unique(relabeled_old_ad)' == 0:(length(ad_nums) - 1)), 'Error in old ad relabeling')
    assert(sum(connected_ad(:)) == (sum(old_ad(:)) + 1),'Error in connected ad finding: %d, %d ',sum(connected_ad(:)),sum(old_ad(:)))
    
    props = regionprops(relabeled_old_ad,'Area','Centroid');

    %if there is only one set of props, we know there was only one adhesion
    %connected to the newest pixel, otherwise, trigger a more complicated
    %procedure
    if (length(props) == 1)
        %pick out the adhesion number from the first entry in find, then
        %check that all other pixels have the same adhesion number and
        %assign the newest pixel to the old adhesion
        ad_number = ad_zamir(find(relabeled_old_ad == 1,1));
        assert(ad_number > 0, 'Error in old ad filtering: adhesion number less than 1');
        assert(all(ad_number == ad_zamir(relabeled_old_ad >= 1)),'Error in old ad filtering: single adhesion with different numbers');

        ad_zamir(pix_pos) = ad_number;
    else
        %Check if all the connected adhesions are below the minimum size,
        %if they are, merge all the adhesions, otherwise, trigger a more
        %complicated procedure
        meets_min = [props.Area] >= min_size;
        if (all([props.Area] < min_size))
            ad_zamir(connected_ad) = min(ad_zamir(old_ad == 1));
        elseif (sum(meets_min) == 1)
            %There is only one large adhesion touching, find the current
            %adhesion number and add the new pixel
            large_area_ad = ismember(relabeled_old_ad,find([props.Area] >= min_size));

            ad_number = ad_zamir(find(large_area_ad == 1,1));
            assert(ad_number > 0, 'Error in large ad filtering: adhesion number less than 1');
            assert(all(ad_number == ad_zamir(large_area_ad == 1)),'Error in large ad filtering: single adhesion with different numbers');

            ad_zamir(connected_ad) = ad_number;
        else
            %There are two or more large adhesions touching the new pixel,
            %figure out which adhesion is closer in terms of centroid to
            %the new pixel and add the new pixel to that adhesion, if there
            %is a tie, just use the first adhesion in the list
            large_area_ads = ismember(relabeled_old_ad,find([props.Area] >= min_size));
            large_area_nums = unique(relabeled_old_ad(large_area_ads));

            %Also, the new pixel might be connected to small adhesions as
            %well as large adhesions, so capture the position of the small
            %adhesions for later assignment to selected large adhesion
            small_area_ads = ismember(relabeled_old_ad, find([props.Area] < min_size));

            closest_relabeled_info = [0,inf];
            for i = 1:length(large_area_nums)
                dist = sqrt((pix_pos_ind(1) - props(i).Centroid(2))^2 + (pix_pos_ind(2) - props(i).Centroid(1))^2);
                if (dist < closest_relabeled_info(2))
                    closest_relabeled_info = [i,dist];
                end
            end
            assert(closest_relabeled_info(1) ~= 0, 'Error in identifying the closest of the large adhesions')

            ad_number = ad_zamir(find(relabeled_old_ad == closest_relabeled_info(1),1));
            assert(ad_number > 0, 'Error in largest ad filtering: adhesion number less than 1');
            assert(all(ad_number == ad_zamir(relabeled_old_ad == closest_relabeled_info(1))),'Error in largest ad filtering: single adhesion with different numbers');

            ad_zamir(pix_pos) = ad_number;
            ad_zamir(small_area_ads) = ad_number;
        end

    end
end

assert(initial_size + 1 == sum(sum(im2bw(ad_zamir,0))), 'Error in adding single pixel: Adhesion set did not grow')
