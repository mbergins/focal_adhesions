cell_mask = imread('../../../../results/focal_adhesions/time_series_01/individual_pictures/001/cell_mask.png');
original_size = size(cell_mask); 
% imwrite(cell_mask,'this_mask.bmp')
% system('/sw/bin/potrace -s this_mask.bmp')
% system('rm this_mask.bmp')
% system('/sw/bin/convert -density 200x200 this_mask.svg this_mask.png')

cell_mask = imread('this_mask.png');
cell_mask = logical(cell_mask(3:end,:,1));

addpath('../../../../src/visualize_cell_features/')

scale_factor = size(cell_mask,1)/original_size(1);

% cell_mask = imresize(cell_mask,scale_factor);

dists = bwdist(not(cell_mask));

green_color_seq = [0.433691756 1.000000000 0.792114695 0.516129032 0.293906810 0.197132616 ...
    0.111111111 0.125448029 0.071684588 0.050179211 0.053763441 0.035842294 ...
    0.043010753 0.035842294 0.057347670 0.039426523 0.046594982 0.043010753 ...
    0.050179211 0.014336918 0.043010753 0.021505376 0.025089606 0.014336918 ...
    0.003584229 0.017921147 0.017921147 0.010752688 0.010752688 0.007168459 ...
    0.000000000 0.000000000 0.003584229 0.000000000 0.003584229];

green_dist_seq = (0:5:5*(length(green_color_seq) - 1))*scale_factor;

green_layer_interp = zeros(size(cell_mask));
for i = 1:size(green_layer_interp,1)
    for j = 1:size(green_layer_interp,2)
        if (dists(i,j) == 0), continue; end
        
        first_color = find(green_dist_seq <= dists(i,j),1,'last');
        if (isempty(first_color)), continue; end
        
        second_color = find(green_dist_seq > dists(i,j),1,'first');
        if (isempty(second_color)), continue; end
        
        assert(second_color == first_color + 1, '%d %d', first_color, second_color);
        
        green_layer_interp(i,j) = interp1([green_dist_seq(first_color),green_dist_seq(second_color)], ...
                                          [green_color_seq(first_color),green_color_seq(second_color)], ...
                                          dists(i,j));
    end
end
disp('Done with Green')

red_color_seq = [0.286738351 0.677419355 0.394265233 0.354838710 0.358422939 0.390681004 ...
    0.430107527 0.365591398 0.318996416 0.293906810 0.207885305 0.182795699 ...
    0.161290323 0.103942652 0.129032258 0.078853047 0.075268817 0.078853047 ...
    0.064516129 0.039426523 0.046594982 0.032258065 0.017921147 0.043010753 ...
    0.007168459 0.010752688 0.014336918 0.017921147 0.014336918 0.021505376 ...
    0.003584229 0.000000000 0.003584229 0.003584229 0.003584229];

red_dist_seq = (0:5:5*(length(red_color_seq) - 1))*scale_factor;

red_layer_interp = zeros(size(cell_mask));
for i = 1:size(red_layer_interp,1)
    for j = 1:size(red_layer_interp,2)
        if (dists(i,j) == 0), continue; end
        
        first_color = find(red_dist_seq <= dists(i,j),1,'last');
        if (isempty(first_color)), continue; end
        
        second_color = find(red_dist_seq > dists(i,j),1,'first');
        if (isempty(second_color)), continue; end
        
        assert(second_color == first_color + 1, '%d %d', first_color, second_color);
        
        red_layer_interp(i,j) = interp1([red_dist_seq(first_color),red_dist_seq(second_color)], ...
                                          [red_color_seq(first_color),red_color_seq(second_color)], ...
                                          dists(i,j));
    end
end
disp('Done with Red')

grad_image = zeros([size(cell_mask),3]);
grad_image(:,:,1) = red_layer_interp;
grad_image(:,:,2) = green_layer_interp;

bounds = find_binary_bounding_box(cell_mask);

for i=1:3
    temp = grad_image(:,:,i);
    
    dilate_mask = imdilate(cell_mask,strel('disk',3));
    
    temp(find(bwperim(dilate_mask))) = 0;
    temp(find(not(dilate_mask))) = 1;
    grad_image(:,:,i) = temp;
end

grad_image = grad_image(bounds(2):bounds(4), bounds(1):bounds(3),:);
imwrite(grad_image,'gradient_image.png');