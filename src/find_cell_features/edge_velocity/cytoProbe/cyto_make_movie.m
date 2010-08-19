

img_path = '/fsm/scripps/analysis/machacek/FocalAdhesions/80505/set_2/pax2_cut1_bleach_corrected/';

filelist_images_b = dir([img_path filesep '*.tif']);
iEntry = 1;
N = length(filelist_images_b);

avimovie = avifile('/home/machacek/movie.avi');

for i=1:N
    if(~filelist_images_b(i).isdir)
        filelist_images(iEntry,:) = filelist_images_b(i).name;
        iEntry = iEntry + 1;
    end
end


cmap = colormap(hot(256));
for i=1:N
    img = imread([img_path filesep filelist_images(i,:)]);
    % scale
    img = double(img);

    img = img ./ 2^14 .* 2^8;
    
    % subtract offset
    img = img-10;
    
    %stretch
    img = img .* 20;
    
    img = uint8(img);
    %M(i) = im2frame(img,cmap);
    f = im2frame(img,cmap);
    avimovie = addframe(avimovie,f);
end
avimovie = close(avimovie);
disp('end');
%movie(img);

%movie2avi(mov,filename'compression'Cinepak')