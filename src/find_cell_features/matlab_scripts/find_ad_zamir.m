function ad_zamir = find_ad_zamir(high_passed_image,binary_image,min_pixel_size,varargin)
% FIND_AD_ZAMIR    Assigns adhesion pixels to specific adhesions using the
%                  same algorithm as described in Zamir, 1999
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Setup variables and parse command line
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

i_p = inputParser;
i_p.FunctionName = 'FIND_AD_ZAMIR';

i_p.addRequired('high_passed_image',@isnumeric);
i_p.addRequired('binary_image',@islogical);
i_p.addRequired('min_pixel_size',@(x)x >= 1);
i_p.addParamValue('debug',0,@(x)x == 1 || x == 0);

i_p.parse(high_passed_image,binary_image,min_pixel_size,varargin{:});

if (i_p.Results.debug == 1), profile on; end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Main Program
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ad_zamir = zeros(size(high_passed_image));

pix_vals = high_passed_image(binary_image);
sorted_pix_vals = sort(unique(pix_vals),'descend');

count = 0;
total_pixels = sum(sum(binary_image));

%Cycle through all pixels of image
for i = 1:length(sorted_pix_vals)
    lin_ind = find(high_passed_image == sorted_pix_vals(i));
    
    
    for j = 1:length(lin_ind)
        [lin_row, lin_col] = ind2sub(size(high_passed_image), lin_ind(j));
        if (not(binary_image(lin_row, lin_col)))
            continue;
        end
    
        assert(ad_zamir(lin_ind(j)) == 0, 'Error: Adhesion already assigned in this position %d',lin_ind(j))
        ad_zamir = add_single_pixel(ad_zamir,lin_ind(j),i_p.Results.min_pixel_size);
        count = count + 1;
    end

    if (mod(count,100) == 0 && i_p.Results.debug)
        disp(['Count: ',num2str(count),'/',num2str(total_pixels)])
    end
end

%renumber the found adhesions to start at one
ad_nums = unique(ad_zamir);
assert(ad_nums(1) == 0, 'Background pixels not found after building adhesion label matrix')
for i = 2:length(ad_nums)
    ad_zamir(ad_zamir == ad_nums(i)) = i - 1;
end

profile off;
if (i_p.Results.debug), profile viewer; end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
%outside the accepted 1-size range
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

%also build an image with the connected old adhesions renumbered to start
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

%if there aren't any pixels which were connected to newest pixel, add the
%newest pixel as a new adhesion, otherwise, start a more complicated
%procedure
if (sum(old_ad(:)) == 0)
    ad_zamir(connected_ad) = max(ad_zamir(:)) + 1;
else
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
