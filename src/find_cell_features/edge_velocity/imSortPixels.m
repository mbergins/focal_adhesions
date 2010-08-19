function sortedPixels = imSortPixels(edge_pixel)
%
%   sorts a chain of pixels 
%
%
%
%
%
%
%
%  Matthias Machacek October 20 2004


% save the initial pixel
init_pix = edge_pixel(1,:);
edge_pixel(1,:) = [];
edge_pixel_sorted = init_pix;

i=1;
i_s=1;
found =1;
while size(edge_pixel,1) > 0 && found == 1
    % Find first 4 connected neighbour
    while i <= size(edge_pixel,1) && sqrt((round(edge_pixel(i,1))-round(edge_pixel_sorted(i_s,1)))^2 +(round(edge_pixel(i,2))-round(edge_pixel_sorted(i_s,2)))^2) > 1
        i=i+1;
    end
    if i>size(edge_pixel,1);
        % Now try to find a 8 connected neighbour
        i=1;
        while  i <= size(edge_pixel,1) && sqrt((round(edge_pixel(i,1))-round(edge_pixel_sorted(i_s,1)))^2 + (round(edge_pixel(i,2))-round(edge_pixel_sorted(i_s,2)))^2) > sqrt(2)
            i=i+1;
        end
        if i>size(edge_pixel,1);
            % No next pixel was found
            found=0;
        end
    end
    if found~=0;
        i_s=i_s+1;
        edge_pixel_sorted(i_s,:) = edge_pixel(i,:);
        edge_pixel(i,:)=[];
        i=1;
    end
end

% check if there is a part missing
if size(edge_pixel,1) > 0
    edge_pixel_sorted2 = init_pix;
    
    i=1;
    i_s=1;
    found =1;
    while size(edge_pixel,1) > 0 && found ==1
        % Find first 4 connected neighbour
        while i <= size(edge_pixel,1) && sqrt((round(edge_pixel(i,1))-round(edge_pixel_sorted2(i_s,1)))^2 +(round(edge_pixel(i,2))-round(edge_pixel_sorted2(i_s,2)))^2) > 1
            i=i+1;
        end
        if i>size(edge_pixel,1);
            % Now try to find a 8 connected neighbour
            i=1;
            while  i <= size(edge_pixel,1) &&  sqrt((round(edge_pixel(i,1))-round(edge_pixel_sorted2(i_s,1)))^2 +(round(edge_pixel(i,2))-round(edge_pixel_sorted2(i_s,2)))^2) > sqrt(2)
                i=i+1;
            end
            if i>size(edge_pixel,1);
                %no next pixel was found
                found=0;
            end
        end
        if found~=0;
            i_s=i_s+1;
            edge_pixel_sorted2(i_s,:) = edge_pixel(i,:);
            edge_pixel(i,:)=[];
            i=1;
        end
    end
    edge_pixel_sorted2(1,:)=[];
    
    
    % connect the two parts
    edge_pixel_sorted = flipud(edge_pixel_sorted);
    l_eps = size(edge_pixel_sorted,1);
    for i=1:size(edge_pixel_sorted2,1)
        edge_pixel_sorted(l_eps+i,:) = edge_pixel_sorted2(i,:);
    end
end
    
sortedPixels = edge_pixel_sorted;


