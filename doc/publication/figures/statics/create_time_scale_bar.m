c_map = jet(198);

bar_width = floor(7500/10);
bar_height = floor(2919/10);

bar = ones(bar_height, bar_width, 3);

alpha = double(ones(bar_height, bar_width));

for i=1:bar_width
    this_height = floor(bar_height*(i/(bar_width+1)));
    
    bar((bar_height - this_height):bar_height,i,1) = c_map(ceil(size(c_map,1)*(i/bar_width)),1);
    bar((bar_height - this_height):bar_height,i,2) = c_map(ceil(size(c_map,1)*(i/bar_width)),2);
    bar((bar_height - this_height):bar_height,i,3) = c_map(ceil(size(c_map,1)*(i/bar_width)),3);
    alpha(1:(bar_height - this_height),i) = 0;
end

imwrite(bar, 'time_scale.png','Alpha',alpha)