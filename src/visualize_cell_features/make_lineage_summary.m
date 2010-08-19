exp_name = 'time_series_01';
tracking_mat = load(fullfile('../../results/focal_adhesions',exp_name,'tracking_matrices/tracking_seq.csv'))+1;
data = load(fullfile('../../results/focal_adhesions',exp_name,'adhesion_props/lin_time_series/Average_adhesion_signal.csv'));
assert(all(size(tracking_mat) == size(data)),'Error: tracking and data matrices are not the same size')

filtered_tracking_mat = zeros(1,size(tracking_mat,2));
filtered_data = zeros(1,size(tracking_mat,2));
for i = 1:size(tracking_mat,1)
    tracking_row = tracking_mat(i,:);
    data_row = data(i,:);
    if (sum(tracking_row >= 1) >= 20)
        filtered_tracking_mat(size(filtered_tracking_mat,1) + 1,:) = tracking_row;
        filtered_data(size(filtered_data,1) + 1,:) = data_row;
    end
end
filtered_tracking_mat = filtered_tracking_mat(2:end,:);
filtered_data = filtered_data(2:end,:);
assert(all(size(filtered_tracking_mat) == size(filtered_data)),'Error: filtered tracking and data matricies are not the same size');

max_data = max(max(filtered_data));
min_data = min(min(filtered_data));

disp(sprintf('Min/Max values %f, %f',min_data,max_data))

cmap = gray(1000);
tracking_mat = filtered_tracking_mat;
image = ones([size(tracking_mat),3]);
for i = 1:size(tracking_mat,1)
    for j = 1:size(tracking_mat,2)
        if (tracking_mat(i,j) <= 0), continue; end
        scaled_data_value = (filtered_data(i,j) - min_data)/(max_data - min_data);
        if (scaled_data_value == 0), scaled_data_value = 1e-10; end
        assert(scaled_data_value <= 1 && scaled_data_value > 0)
        
        cmap_row = ceil(scaled_data_value*size(cmap,1));
        assert(cmap_row > 0 && cmap_row <= size(cmap,1),'Cmap row is out of range: %d', cmap_row);
        
        image(i,j,1:3) = cat(3,cmap(cmap_row,1),cmap(cmap_row,2),cmap(cmap_row,3));
    end
end

scale = zeros(10,size(image,2),3);
for i = 1:size(image,2)
    cmap_row = ceil((i/size(image,2))*size(cmap,1));
    assert(cmap_row > 0 && cmap_row <= size(cmap,1),'Cmap row is out of range: %d', cmap_row);
    
    scale(:,i,1:3) = repmat(cat(3,cmap(cmap_row,1),cmap(cmap_row,2),cmap(cmap_row,3)), [size(scale,1), 1]);
end

composite_image = [scale;repmat(1,[1,size(image,2),3]);image];


imwrite(composite_image,['lineage_sum_',exp_name,'.png'])