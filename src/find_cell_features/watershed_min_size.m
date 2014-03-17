function final_label_mat = watershed_min_size(imageOrig,threshed_image,min_size)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Faster version of watershed min size algorithm from Hoffman Lab
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

im = imageOrig.*(threshed_image > 0);      %threshold (should eliminate any negative values)

%reset any values at or below zero to a small positive value to ensure that
%the find further down includes all segmented pixels
im(threshed_image > 0 & im <= 0) = 1E10;   

%The following code might try to lookup a pixel value outside the image
%size, so add a layer of zeros to deal with that possibility. We will
%remove the layer at the end of processing.
mat = padarray(im, [1 1]);     %pads matrix with 0s on all sides

% im_label = bwlabel(im > 0,8);
% props = regionprops(im_label,'Area');
% small_FA = ismember(im_label,find([props.Area] < min_size*2));
% labelMat = bwlabel(small_FA,8);
% im(small_FA) = 0;

[r,c,v] = find(mat);          %collect non-zero values [v] and their location [r,c]
list = [r c v];               %pacakge output vectors into one matrix
list = sortrows(list, -3);    %sorts rows by v (brightest to dimmest)

%% Identify Patches

labelMat = zeros(size(mat));  %pre-allocate matrix to collect patch numbers
patchNum = max(labelMat(:)) + 1;   %independent index for patch labels
for i = 1:size(list,1)
    hood = getFour(list(i,1:2), labelMat);  %look in 'hood for existing patches
    patchList = unique(nonzeros(hood));      %find unique, non-zero patch labels
    
    switch length(patchList)
        case 0                             %no patches in 'hood
            labelMat(list(i,1),list(i,2)) = patchNum; %assign new patch number
            patchNum = patchNum+1;
        case 1                             %one neighboring patch
            labelMat(list(i,1),list(i,2)) = patchList; %assign to the existing patch
        otherwise                          %>1 neighboring patch\
            %We need to know which touching patches are over the minimum
            %size criteria and while we are looking at them, we will also
            %gather patch intensity information.
            allInd = []; aver_int = []; sz = [];
            for j = 1:length(patchList)    %for each patch in the list
                ind = find(labelMat == patchList(j)); %find all pixels with corresponding patch number
                sz(j) = length(ind);       %#ok<AGROW> %patch size
                aver_int(j) = sum(sum(mat(ind)))/sz(j);    %#ok<AGROW> %patch integrated intensity
                allInd = [allInd; ind];    %#ok<AGROW> %collect all indicies
            end
            
            %This bit of code finds the index in the patchList that has the
            %highest intensity, but only considers the patch numbers with
            %enough size to avoid being merged. This patch number is then
            %merged/assigned to the current pixel or other adjacent
            %patches.
            [~, brightest_large_patch_index] = max((sz >= min_size) .* aver_int);
            brightest_large_patch_num = patchList(brightest_large_patch_index);
            
            doesnt_meet_size_nums = patchList(sz < min_size);
            for small_patch_num = doesnt_meet_size_nums'
                labelMat(labelMat == small_patch_num) = brightest_large_patch_num;
            end
            labelMat(list(i,1),list(i,2)) = brightest_large_patch_num;
    end
end

final_label_mat = labelMat(2:(size(labelMat,1)-1), 2:(size(labelMat,2)-1)); %remove padding

%Renumber the adhesions using the original threshed image as a guide, this
%will ensure that the numbers increase left to right and up to down.
no_split_labeling = bwlabel(threshed_image,4);
ad_renumber = zeros(size(final_label_mat));
new_ad_nums = 1;
for ad_num = 1:max(no_split_labeling(:))
    split_ad_nums = unique(final_label_mat(no_split_labeling == ad_num));
    
    for j=1:length(split_ad_nums)
        if (length(split_ad_nums) > 1)
            1;
        end
        ad_renumber(final_label_mat == split_ad_nums(j)) = new_ad_nums;
        new_ad_nums = new_ad_nums + 1;
    end
end

final_label_mat = ad_renumber;
return

%% Sub Functions

function hood = getEight(index, mat)
% finds pixels with eight-point connectivity to the input pixel
r = index(1);
c = index(2);

hood = [
    mat(r-1,c-1)
    mat(r-1,c  )
    mat(r-1,c+1)
    mat(r  ,c-1)
    mat(r  ,c+1)
    mat(r+1,c-1)
    mat(r+1,c  )
    mat(r+1,c+1)
    ];
return

function hood = getFour(index, mat)
% finds pixels with four-point connectivity to the input pixel
r = index(1);
c = index(2);

hood = [
    mat(r-1,c  )
    mat(r  ,c-1)
    mat(r  ,c+1)
    mat(r+1,c  )
    ];
return

