exp_name = 'time_series_01';
tracking_mat = load(fullfile('../../results/focal_adhesions',exp_name,'tracking_matrices/tracking_seq.csv'))+1;
data = load(fullfile('../../results/focal_adhesions',exp_name,'adhesion_props/lin_time_series/Angle_to_center.csv'));
% orient_data = load(fullfile('../../results/focal_adhesions',exp_name,'adhesion_props/lin_time_series/Orientation.csv'));
area_data = load(fullfile('../../results/focal_adhesions',exp_name,'adhesion_props/lin_time_series/Area.csv'));
assert(all(size(tracking_mat) == size(data)),'Error: tracking and data matrices are not the same size')

% ad_lifetimes = zeros(size(tracking_mat,1),1);
% for i = 1:size(tracking_mat,1), ad_lifetimes(i) = sum(tracking_mat(i,:) >= 1); end
% 
% ad_lifetime_remain = zeros(size(tracking_mat));
% for i = 1:size(tracking_mat,1)
%     tracking_row = tracking_mat(i,:);
%     for j = 1:size(tracking_mat,2)
%         ad_lifetime_remain(i,j) = sum(tracking_row(j:end) >= 1);
%     end
% end

%Build the area sums matrix, summing up all adhesions at a given angle from
%the cell centroid, based on breaking the angles into as many discrete
%region as the number of rows in the 'area_sums' variable
area_sums = zeros(360*2,size(tracking_mat,2));
for i = 1:size(tracking_mat,1)
    tracking_row = tracking_mat(i,:);
    for j = 1:size(tracking_mat,2)
        if (isnan(data(i,j))), continue; end

        angle = data(i,j) * (360/(2*pi));
        if (angle == 0), angle = 360; end
        assert(angle <= 360 & angle > 0,'%d %d',angle, data(i,j));

        scaled_angle = ceil(size(area_sums,1) * (angle/360));
        assert(scaled_angle > 0 & scaled_angle <= size(area_sums,1));

        area_sums(scaled_angle,j) = area_sums(scaled_angle,j) + area_data(i,j);
    end
end
area_sums = log(area_sums);
disp(sprintf('Min/Max Log Area Sums: %f %f',min(min(area_sums(area_sums ~= -Inf))), max(area_sums(:))));

image = zeros(size(area_sums,1), size(area_sums,2), 3);
cmap = hot(1000);

min_area_sum = min(min(area_sums(area_sums ~= -Inf)));
range_area_sum = range(area_sums(area_sums ~= -Inf));

for i = 1:size(area_sums,1)
    for j = 1:size(image,2)
        if (area_sums(i,j) == -Inf), continue; end
        
        scaled_area = (area_sums(i,j) - min_area_sum)/range_area_sums;
        if (scaled_area == 0), scaled_area = 1e-10; end
        assert(scaled_area > 0 && scaled_area <= 1, '%d %d', scaled_area, area_sums(i,j));
        
        cmap_row = ceil(scaled_area*size(cmap,1));
        assert(cmap_row > 0 && cmap_row <= size(cmap,1),'Cmap row is out of range: %d', cmap_row);
        
        image(i,j,1:3) = cat(3,cmap(cmap_row,1),cmap(cmap_row,2),cmap(cmap_row,3));
    end
    if (mod(i,10) == 0), disp(num2str(i)); end
end

scale = zeros(10,size(image,2),3);
for i = 1:size(image,2)
    cmap_row = ceil((i/size(image,2))*size(cmap,1));
    assert(cmap_row > 0 && cmap_row <= size(cmap,1),'Cmap row is out of range: %d', cmap_row);
    
    scale(:,i,1:3) = repmat(cat(3,cmap(cmap_row,1),cmap(cmap_row,2),cmap(cmap_row,3)), [size(scale,1), 1]);
end

composite_image = [image;repmat(0.5,[1,size(image,2),3]);scale];

imwrite(composite_image,['as_',exp_name,'.png'])