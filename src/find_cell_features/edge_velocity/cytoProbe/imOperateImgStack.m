function imOperateImgStack(data, operation)
%
%
% get GUI handle
hfsmC=findall(0,'Tag','cytoProbeFigure');
handlesCytoProbe = guidata(hfsmC);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%       Load the data       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if operation == 1 || operation == 2 || operation == 7 ||...
   operation == 13
    % 1: multiplication
    % 2: division
    % 7: correlation
    % 13: bleach correction
    [fileList1]=getFileStackNamesAdv(data.image1_path,data.read_mode);
    [fileList2]=getFileStackNamesAdv(data.image2_path,data.read_mode);
elseif  operation == 3 
    % 3:  
    [fileImageList1]=getFileStackNamesAdv(data.image1_path,data.read_mode); 
    % load av_intensity
    load(data.sVector_path);
elseif  operation == 4 || operation == 8 || operation == 11 || ...
                          operation == 18 || operation == 19
    % 4:
    % 8: 
    % 11: 
    % 18: generate avi movie
    % 19: calculate average image  
    [fileList1]=getFileStackNamesAdv(data.image1_path,data.read_mode);  
elseif operation == 5  
    % 5: FRET  
    [fileList1]=getFileStackNamesAdv(data.image1_path,data.read_mode);
    [fileList2]=getFileStackNamesAdv(data.image2_path,data.read_mode); 
    if ~strcmp(data.image3_path,'-- Choose an scaling vector --')
        % get the YFP channel
        [fileList3]=getFileStackNamesAdv(data.image3_path,data.read_mode);
        YFP_correction = 1;
    else
        YFP_correction = 0;
    end
elseif  operation == 6 || operation == 26
    % 6  :  get activity from edge distance
    % 26 :  Do PKA FRET (UCSD)
    [fileList1]=getFileStackNamesAdv(data.image1_path,data.read_mode);  
    [fileList2]=getFileStackNamesAdv(data.image2_path,data.read_mode);
elseif  operation == 9  
    % subtract images
    [d1, d2, ext1] = fileparts(data.image1_path);
    [d1, d2, ext2] = fileparts(data.image2_path);    
    if strcmp(ext2,'.mat') & ~strcmp(ext1,'.mat')
        [fileList1]=getFileStackNamesAdv(data.image1_path,data.read_mode);
        load(data.image2_path);
        vector_mode = 1;  
    else
        [fileList1]=getFileStackNamesAdv(data.image1_path,data.read_mode);
        [fileList2]=getFileStackNamesAdv(data.image2_path,data.read_mode);
        vector_mode = 0;
    end
elseif  operation == 10
    %subtract vector field
    [fileList1]=getFileStackNamesAdv(data.image1_path,data.read_mode);
    [fileList2]=getFileStackNamesAdv(data.image2_path,data.read_mode);
elseif  operation == 14
    % get curvature from spline
    splines = load(data.image1_path); 
    a= fieldnames(splines);
    a1 = char(a{1});
    a2 = char(a{2}); 
    x_spline = splines.(a1);
    y_spline = splines.(a2);
elseif  operation == 15
    % subtract bg
    [fileList1]=getFileStackNamesAdv(data.image1_path,data.read_mode);
    if ~strcmp(data.image2_path,'-- Choose an image --')
        [fileList2]=getFileStackNamesAdv(data.image2_path,data.read_mode);   
        roi_bg = 0;
    else
        % user must define RIO for bg
        roi_bg = 1;
    end
elseif  operation == 16
    % translate image
    [fileList1]=getFileStackNamesAdv(data.image1_path,data.read_mode);  
elseif  operation == 17
    % correct for transformation
    [fileList1]=getFileStackNamesAdv(data.image1_path,data.read_mode);    
    [fileList2]=getFileStackNamesAdv(data.image2_path,data.read_mode); 
elseif  operation == 20
    % overlay pixels
    [fileList1]=getFileStackNamesAdv(data.image1_path,data.read_mode);    
    pixels_tmp = load(data.image2_path); 
    a = char(fieldnames(pixels_tmp));
    pixels = pixels_tmp.(a);
elseif  operation == 21    
    % average  vector field
    [fileList1]=getFileStackNamesAdv(data.image1_path,data.read_mode);
%     if ~strcmp(data.image2_path,'-- Choose an image --')
%         [fileList2]=getFileStackNamesAdv(data.image2_path,data.read_mode);   
%         roi_bg = 0;
%     else
%         % user must define RIO for bg
%         roi_bg = 1;
%     end
elseif  operation == 22 | 23 | 24| 25
    % 22: convert MPM matrix to a of vector files
    % 23: analyse speckle track life time and speeds
    % 24: convert SCORES matrix to a stack of vector files   (kinScores
    % format
    % 25: analyse protrusion  
    [fileList1]=getFileStackNamesAdv(data.image1_path,data.read_mode);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%   End load the data       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Create a directory for the results
if operation == 1
    resultDir = [data.result 'multiplied'];
    file_app = 'multiplied_';
elseif  operation == 2
    resultDir = [data.result 'divided'];    
    file_app = 'divided_';
elseif   operation == 3
    resultDir = [data.result 'scaled'];    
    file_app = 'scaled_';
elseif operation == 4
    resultDir = data.result;  
elseif operation == 5
    resultDir = [data.result 'FRETn'];  
    file_app = 'fret_c_n_';
elseif operation == 6
    resultDir = [data.result 'edge_analysis'];  
    file_app = 'activity_';  
elseif operation == 7
    resultDir = [data.result 'correlations'];  
    file_app = 'corr_';     
elseif operation == 8
    resultDir = [data.result 'inverse'];  
    file_app = 'inv_';  
elseif operation == 9
    resultDir = [data.result 'subtract'];  
    file_app = 'difference_';     
elseif operation == 10
    resultDir = [data.result 'subtract_v'];  
    file_app = 'difference_v_';    
elseif operation == 11
    resultDir = [data.result 'cell_mask'];  
    file_app = 'mask_';   
elseif operation == 13
    resultDir = [data.result 'bleach_corrected'];  
    file_app = 'bg_c_';    
elseif operation == 14
    resultDir = [data.result 'edge_curvature'];  
    file_app = 'curvature_';   
elseif operation == 15
    resultDir = [data.result 'bg_subtracted'];  
    file_app = 'bg_c_';  
elseif operation == 16
    resultDir = [data.result 'translated'];  
    file_app = 't_';  
elseif operation == 17
    resultDir = [data.result 'translation_corrected'];  
    file_app = 'tc_';    
elseif operation == 18
    resultDir = [data.result 'movie'];  
    file_app = 'mov_';    
elseif operation == 19
    resultDir = [data.result 'average_image'];  
    file_app = 'av_';   
elseif operation == 20
    resultDir = [data.result 'overlay_pixels'];  
    file_app = 'ov_';   
elseif operation == 21
    resultDir = [data.result 'average_vector_field'];  
    file_app = 'av_vec_field_'; 
elseif operation == 22
    resultDir = [data.result 'MPM_vectors'];  
    file_app = 'vector_';   
elseif operation == 23
    resultDir = [data.result 'track_analysis'];  
    file_app = 'vector_'; 
elseif operation == 24
    resultDir = [data.result 'SCORES_vectors'];  
    file_app = 'scorea_';   
elseif operation == 25
    resultDir = [data.result 'protrusion_analysis'];  
    file_app = 'scorea_';   
elseif operation == 26
    resultDir = [data.result 'YFP_CFP_ratio'];  
    file_app = 'YFP_CFP_ratio_';       
end

[s, mess, messid] = mkdir(resultDir);
if s==0
    disp('Failed to create the masked images directory');
    return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%  Get the number of files to process 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
if operation == 1 
    options         = cytoMultiplicationDlg;    
    % use only first image 2 to multiply entire image 1 stack 
    first_image = options.first_image;
    % stretch image 
    stretch = options.stretch;
    stretch_i = options.stretch_i;
    [d nr_images1] = size(fileList1);
    [d nr_images2] = size(fileList2);
    if first_image
        nr = nr_images1;
    else
        nr = min(nr_images1, nr_images2);   
    end    
elseif operation == 2
    options         = cytoDivisionDlg;    
    % use only first image 2 to divide entire image 1 stack 
    first_image = options.first_image;
    % stretch image 
    stretch = options.stretch;
    stretch_i = options.stretch_i;
    stretch_c = options.stretch_c;
    
    if stretch_c 
        stretch_coeff = options.stretch_coeff;
    end
    [d nr_images1] = size(fileList1);
    [d nr_images2] = size(fileList2);
    if first_image
        nr = nr_images1;
    else
        nr = min(nr_images1, nr_images2);   
    end
elseif operation == 3
    [d nr_images1] = size(fileList1);
    nr_sVector = length(av_intensity_norm_fit);
    
    nr = min(nr_images1, nr_sVector);
elseif operation == 6 || operation == 7 || operation == 13
    [d nr_images1] = size(fileList1);
    [d nr_images2] = size(fileList2);

    nr = min(nr_images1, nr_images2);   
elseif operation == 5 
    [d nr_images1] = size(fileList1);
    [d nr_images2] = size(fileList2);
    if YFP_correction
        [d nr_images3] = size(fileList3);
        nr = min([nr_images1, nr_images2, nr_images3]); 
    else
        nr = min(nr_images1, nr_images2); 
    end
elseif operation == 9   
    if vector_mode == 0
        [d nr_elements1] = size(fileList1);
        [d nr_elements2] = size(fileList2);
        if nr_elements2 == 1
            % switch to one image mode
            nr = nr_elements1;
            one_image_mode = 1;
        else
            nr = min(nr_elements1, nr_elements2);   
            one_image_mode = 0;
        end
    elseif vector_mode == 1
        [d nr_elements1] = size(fileList1);
        nr_sVector = length(av_intensity);
        nr = min(nr_images1, nr_sVector);    
    end
elseif operation == 10         
    [d nr_elements1] = size(fileList1);
    [d nr_elements2] = size(fileList2);
    nr = min(nr_elements1, nr_elements2);   
elseif operation == 11         
    [d nr_elements1] = size(fileList1);
    nr = nr_elements1;     
elseif operation == 14
     nr = size(x_spline,2);
elseif operation == 15
     [d nr_elements1] = size(fileList1);
     nr = nr_elements1;
elseif operation == 16
     [d nr_elements1] = size(fileList1);
     nr = nr_elements1;     
elseif operation == 17
    [d nr_elements1] = size(fileList1);
    [d nr_elements2] = size(fileList2);
    nr = min(nr_elements1, nr_elements2);  
elseif operation == 18
    [d nr_elements1] = size(fileList1);
    nr = nr_elements1;  
elseif operation == 19
    [d nr_elements1] = size(fileList1);
    nr = nr_elements1;  
elseif operation == 20
    [d nr_elements1] = size(fileList1);
    nr = nr_elements1;      
elseif operation == 21
    [d nr_elements1] = size(fileList1);
    nr = nr_elements1;     
elseif operation == 22
    nr = 1;     
elseif operation == 23
    nr = 1;     
elseif operation == 24
    nr = 1;  
elseif operation == 25
    nr = 1;    
elseif operation == 26
    [d nr_elements1] = size(fileList1);
    nr = nr_elements1;   
else
    [d nr] = size(fileList1); 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%   Get some specific settings
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if operation == 5
    % Rac FRET
    options           = cytoRacFRETDlg(0.395, 0.06); %(handles.CFP_bleedthrough, handles.YFP_bleedthrough);    
    background        = options.background;    
    CFP_bleedthrough  = options.CFP_bleedthrough;
    YFP_bleedthrough  = options.YFP_bleedthrough;
    if YFP_bleedthrough == 0
        YFP_correction = 0
    end
    %handles.CFP_bleedthrough  = CFP_bleedthrough;
    %handles.YFP_bleedthrough  = YFP_bleedthrough;
    
    if background
        img = imread(fileList1{1});    
        img = imadjust(img,stretchlim(img),[]);
        set(handlesCytoProbe.info,'String',['Choose representative background region. Finish with double click', resultDir]);
        %cytoROIinfo;
        figure
        roi_mask = roipoly(img);        
        % Get average background
        for i=1:nr
            img1 = imread(fileList1{i});  
            img2 = imread(fileList2{i});  
            bg_img1 = double(img1) .* double(roi_mask);
            bg_img2 = double(img2) .* double(roi_mask);
            bg1(i) = mean(bg_img1(find(bg_img1>0)));
            bg2(i) = mean(bg_img2(find(bg_img2>0)));
            if YFP_correction
                img3 = imread(fileList3{i});
                bg_img3 = double(img3) .* double(roi_mask);
                bg3(i) = mean(bg_img3(find(bg_img3>0)));
            end
        end
        av_bg1 = mean(bg1);
        av_bg2 = mean(bg2);
        if YFP_correction
            av_bg3 = mean(bg3);
        end
    end
end

if operation == 6
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%   Activity from edge   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    options = cytoBandDlg;
    band_width      = options.band_width;
    background      = options.background;
    roi             = options.roi;
    scores_data     = options.scores;
    vectors_data    = options.vectors;
    images_data     = options.images;
    nr_user_bands   = options.nr_user_bands;
    normalise       = options.normalise;
    scores_values   = options.scores_values;
    % :0 use all scores values
    % :1 use only positive scores values
    % :2 use only negative scores values
    if  band_width <= 0
        return;
    end
    if roi
        if images_data
            img = imread(fileList1{1});    
            img=imadjust(img,stretchlim(img),[]);
        else
            img = imread(fileList2{1});    
            %img = imadjust(img,stretchlim(img),[]);           
            %vecMap1 = load(fileImageList1{i});
        end
        set(handlesCytoProbe.info,'String',['Choose region of interest:  ', resultDir]);
        %cytoROIinfo;
        figure
        roi_mask = roipoly(img);
    end
end

if operation == 15
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%   Subtract background   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%      
    if ~roi_bg
        options = cytoSubBgDlg;
        if options.roi_bg
            set(handlesCytoProbe.info,'String',['Choose background region:  ', resultDir]);
            img = imread(fileList1{1});
            img=imadjust(img,stretchlim(img),[]);
            figure
            roi_mask = roipoly(img);
        end
    else
        set(handlesCytoProbe.info,'String',['Choose background region:  ', resultDir]);
        img = imread(fileList1{1});
        img=imadjust(img,stretchlim(img),[]);
        figure
        roi_mask = roipoly(img);       
    end
end

if operation == 16
    options = cytoTranslationDlg(0,0);
    x_shift = options.x_shift;
    y_shift = options.y_shift;
end

if operation == 18
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%   Make movie    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
    clear aviobj;
    % check for system
    if isunix
        % initialize movie
        avimovie = avifile([resultDir filesep 'movie.avi']);        
    else
        % initialize movie
        avimovie = avifile([resultDir filesep 'movie.avi'],'compression','Cinepak');
    end
    % get options
    options = cytoMovieDlg;
    
    % set frames per second
    avimovie.fps = options.fps;
    
    % check if to draw a box
    if options.box
        h = figure;
        img = imread(fileList1{1});
        img = imadjust(img);
        if options.img_resize ~= 1
            img = imresize(img,options.img_resize,'bicubic');
        end
        imshow(img);
        title('Choose upper left and lower right corners of the box.');
        box_coord = ginput(2);
        box_coord = int16(box_coord);
        close(h);
    end
    
    % check if to draw a bar
    if options.bar
        h = figure;
        img = imread(fileList1{1});
        if ndims(img) < 3
            % adjust the contrast of the image
            img = imadjust(img);
        end
        if options.img_resize ~= 1
            img = imresize(img,options.img_resize,'bicubic');
        end
        imshow(img);
        title('Choose position of the bar.');
        bar_coord = ginput(1);
        bar_coord = int16(bar_coord);
        close(h);        
        bar_length = options.bar_length;
            
        % adjust bar to the image scaling
        bar_length = bar_length .* options.img_resize;
    end    
    
    % Add number
    if options.timer
        % Get position where to add text
        img = imread(fileList1{1});   
        h = figure;
        imshow(img,[]);
        title('Choose position of the text.');
        text_coord = ginput(1);
        text_coord(2)= size(img,1)- text_coord(2);
        text_coord = int16(text_coord);
        
        close(h);        
        
        % Make an image the same size and put text in it 
        for i_img = 1:nr
            hf = figure;
            image(ones(size(img)));
            set(gca,'units','pixels','position',[5 5 size(img,2)-1 size(img,1)-1],'visible','off')
            
            time_string = sprintf('%05d',  options.timer_interval * (i_img - 1));
            text_string = [time_string, ' s'] ;
            % Text at arbitrary position
            text('units','pixels','color','white','position',text_coord,'fontsize',12,'string',text_string)

            % Capture the text image
            % Note that the size will have changed by about 1 pixel
            tim = getframe(gca);
            close(hf)
            
            % make sure tim.cdata has the right size
            tim.cdata(size(img,1):end, :,:) = [];
            tim.cdata(:, size(img,2):end,:) = [];
            tim.cdata(size(img,1), :,:) = 0;
            tim.cdata(:, size(img,2),:) = 0;            
            
            % Extract the cdata
            text_img = tim.cdata(:,:,1);

            % Make a mask with the negative of the text
            [i_text,j_text] = find(text_img > 0);

            % Get the indicies of the text
            ind_text{i_img} = sub2ind(size(img(:,:,1)),i_text,j_text);
            %ind_text{i_img} = sub2ind(size(text_img),i_text,j_text);
        end
    end
end

if operation == 19
    % average images
    img_tmp = imread(fileList1{1}); 
    sum_img = zeros(size(img_tmp));
    counter_img = zeros(size(img_tmp));
end

if operation == 20
    % overlay pixels
end

if operation == 21
    % average vector field in ROI
    % get image to set ROI
    [FileName, PathName] = uigetfile;
    img = imread([PathName, FileName]);
    set(handlesCytoProbe.info,'String',['Choose region of interest:  ', resultDir]);
    %cytoROIinfo;
    h = figure;
    [roi_mask, x_roi, y_roi] = roipoly(img);
    close(h);
end

if operation == 26
    % PKA FRET
    % options           = cytoRacFRETDlg(0.395, 0.06); %(handles.CFP_bleedthrough, handles.YFP_bleedthrough);    
    % UCSD options
    % call interface to determine bleedthrough factors
    options = cytoRacFRETDlg(0.3, 0.01); %(handles.CFP_bleedthrough, handles.YFP_bleedthrough);    
    
    % get flatfield
    flatfield_img = imread(fileList2{1});
 
    % filter flatfield image
    h = fspecial('gaussian',5,0.5);
    flatfield_img =  medfilt2(flatfield_img,[3 3]);
    flatfield_img =  imfilter(flatfield_img,h);
    
    background        = options.background;    
    CFP_bleedthrough  = options.CFP_bleedthrough;
    YFP_bleedthrough  = options.YFP_bleedthrough;

    
    if background
        img = imread(fileList1{1});    
        img = imadjust(img,stretchlim(img),[]);
        set(handlesCytoProbe.info,'String',['Choose representative background in region. Finish with double click', resultDir]);
        %cytoROIinfo;
        figure
        roi_mask = roipoly(img);
        
        figure
        roi_mask = roipoly(img);    
        % Get average background
        for i=1:nr
            img1 = imread(fileList1{i});  

            bg_img1 = double(img1) .* double(roi_mask);
 
            bg1(i) = mean(bg_img1(find(bg_img1>0)));
            bg2(i) = mean(bg_img2(find(bg_img2>0)));
            if YFP_correction
                img3 = imread(fileList3{i});
                bg_img3 = double(img3) .* double(roi_mask);
                bg3(i) = mean(bg_img3(find(bg_img3>0)));
            end
        end
        av_bg1 = mean(bg1);
        av_bg2 = mean(bg2);
        if YFP_correction
            av_bg3 = mean(bg3);
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% allocate memory
if  operation == 1 || operation == 2
    img = imread(fileList1{1});
    [w, h] = size(img);
    %res_image = zeros(w,h,nr); 
    res_image = cell(1,nr);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%          Main loop          %%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
warning off all
intwarning('off')
nr
strg = sprintf('%%.%dd',5);
backSpc = ['\b\b\b\b\b'];
for i=1:nr

    fprintf(1,[strg],i); 
    
    if operation == 1
        % multiply
        img1 = imread(fileList1{i});
        if i > 1 
            if ~first_image
                img2 = imread(fileList2{i}); 
            end
        else
            img2 = imread(fileList2{i});          
        end
        res_image{1,i} = double(img1) .* double(img2);
        max_img_val(i) = max(max(res_image{1,i}));
    elseif  operation == 2
        % divide
        img1 = imread(fileList1{i});
        if i > 1 
            if ~first_image
                img2 = imread(fileList2{i}); 
            end
        else
            img2 = imread(fileList2{i});          
        end
        res_image{1,i} = double(img1) ./ double(img2);
        max_img_val(i) = max(max(isfinite(res_image{1,i}).* res_image{1,i}));
    elseif operation == 3
        img1 = imread(fileList1{i});
        %res_image = immultiply(img1, av_intensity); 
        res_image = uint16(double(img1) .* av_intensity_norm_fit(i)); 
    elseif operation == 4
        % average intensity of non-zero pixels
        img1 = imread(fileList1{i}); 
        mask = img1 > 0;
        nr_pix = sum(sum(mask));
        total_intensity = double(sum(sum(img1)));
        av_intensity(i) = total_intensity / nr_pix;
    elseif operation == 5
        img1 = imread(fileList1{i}); 
        img2 = imread(fileList2{i});
        % if there is a YFP image, correct for YFP bleedthrough
        if YFP_correction
            img3 = imread(fileList3{i});
        end
        % corrected an normalized fret images as used for Rac
        if background 
            img1_c = double(img1) - av_bg1;
            img2_c = double(img2) - av_bg2;
            if YFP_correction
                img3_c = double(img3) - av_bg3;
            end
        else
            img1_c = double(img1);
            img2_c = double(img2);  
            if YFP_correction
                img3_c = double(img3);
            end
        end
        % calculation of the normalized and corrected FRET image
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if YFP_correction
            res_image = (img1_c - CFP_bleedthrough .* img2_c - YFP_bleedthrough .* img3_c) ./ img2_c; 
        else
            res_image = (img1_c - CFP_bleedthrough .* img2_c) ./ img2_c;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % get rid of minus values outside of cell
        res_image = res_image .* (res_image >= 0);
        % get rid of values > 1 outside of cell
        %res_image = res_image .* (res_image < 1);
        %res_image = uint16(res_image .* 2^16);
        res_image = uint16(res_image .* 1000);
        
        
        
        
        
    elseif operation == 6
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Calculate values as function of distance from the cell edge  %%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % generate bands   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % generate the distance map from the mask
        img2 = imread(fileList2{i});
        distance_image = bwdist(~img2);

        % generate an offset off edge_off_set to accound for data
        % outside of the edge
        if 0
            edge_off_set = 31;
            se = strel('disk',edge_off_set);
            img2 = imerode(img2,se);
        end
        
        % generate regions bands
        res_image = zeros(size(distance_image));
        % get the number of bands
        if roi
            distance_image = roi_mask .* distance_image;
        end
        nr_bands = floor(max(distance_image(:)) / band_width) - 3;
        if nr_bands > nr_user_bands & nr_user_bands ~= 0
            nr_bands = nr_user_bands;
        end
        % end generate bands   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        if scores_data
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%          Scores                %%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            s_list = load(fileList1{i});
            a=char(fieldnames(s_list));
                
            %no_zero_index = find(s_list.(a)(:,2));
            %scores = s_list.(a)(no_zero_index,:);
            scores = s_list.(a);
            
            if scores_values == 1
                % take only positive scores
                scores = scores(find(scores(:,4) > 0),:);
            elseif scores_values == 2
                % take only negative scores
                scores = scores(find(scores(:,4) < 0),:);
            end
            
            if 1
                if ~isempty(scores) & scores(1,3) ~= 0
                    for j_b=1:nr_bands
                        band_img = (distance_image >= (j_b-1)*band_width) - (distance_image > j_b*band_width);
                        nr_sc = 0;
                        scores_sum = 0;
                        
                        for v=1:size(scores,1)
                            if band_img(round(scores(v,2)),round(scores(v,3)))
                                scores_sum = scores_sum + scores(v,4);
                                nr_sc = nr_sc + 1;
                            end
                        end
                        band_area = sum(band_img(:));
                        av_activity_band{i}(j_b) = scores_sum / band_area;
                        score_num{i}(j_b) = nr_sc;
                        score_density{i}(j_b) = score_num{i}(j_b)/band_area;
                        %av_activity_band{i}(j_b)=mean(scores(in_shift,4));
                    end
                end
            else            
                if ~isempty(scores)
                    for j_b=1:nr_bands
                        % generate bands
                        band = (distance_image >= (j_b-1)*band_width) - (distance_image > j_b*band_width);
                        if roi
                            band = band .* roi_mask; 
                        end
                        band_polygon = cell2mat(bwboundaries(band, 'noholes'));
                        
                        % delet every second element of band to make code run
                        % faster
                        p_length = size(band_polygon,1);
                        del_ind = 1:3:p_length;
                        band_polygon(del_ind,:) = [];

                        %get the area of the band
                        band_area = polyarea(band_polygon(:,2),band_polygon(:,1));

                        % find scores located inside band
                        in_shift = inpolygon(scores(:,3),scores(:,2),band_polygon(:,2),band_polygon(:,1));
                        av_activity_band{i}(j_b)=sum(scores(find(in_shift),4)) / band_area;

                        score_num{i}(j_b) = length(find(in_shift));
                        score_density{i}(j_b) = score_num{i}(j_b)/band_area;
                    end
                end
            end
            clear scores;
            %figure, plot(score_data_extr(:,1),score_data_extr(:,2),'x');
        elseif images_data
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%         Images                 %%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % read the activity image
            img1 = imread(fileList1{i});
             
            % get average background intensity
            if background
                inv_bg = double(~img2);
                img1_bg = inv_bg .* double(img1);
                av_img1_bg = sum(img1_bg)/sum(inv_bg);
                img1_bg_corr = double(img1)-av_img1_bg;
                img1 = img1_bg_corr.* (img1_bg_corr >= 0);
            end
            
            % extract data in band
            for j_b=1:nr_bands
                band = (distance_image >= (j_b-1)*band_width) - (distance_image > j_b*band_width);
                if roi
                   band = band .* roi_mask; 
                end
                if mod(j_b,2) ~= 0
                    res_image = res_image + bwmorph(band, 'remove');
                end
                % extract data in band
                activity_band = band .* double(img1);
                % average activity in band
                av_activity_band{i}(j_b) = sum(sum(activity_band,1),2) / length(find(activity_band));
            end
            img1_temp = nrm(img1,1); 
            img1_temp(find(res_image)) = 1;
            res_image = img1_temp;
        elseif vectors_data
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%           Vectors              %%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % read the vector map
            vecMap1 = load(fileList1{i});
                   
            for j_b=1:nr_bands
                % extract data in band   
                av_x_v = 0;
                av_y_v = 0;
                l_vec = 0;
                nr_vec = 0;
                a=char(fieldnames(vecMap1));
                vectors = vecMap1.(a);
                
                band_img = (distance_image >= (j_b-1)*band_width) - (distance_image > j_b*band_width);
                         
                % eliminate nan in vector
                vectors(isnan(vectors(:,3)),:) = [];
                
                for v=1:size(vectors,1)      
                    if band_img(round(vectors(v,1)),round(vectors(v,2)))
                            y_v = vectors(v,3) - vectors(v,1);
                            x_v = vectors(v,4) - vectors(v,2);
                            l_vec = l_vec+sqrt(x_v^2 + y_v^2);
                            av_y_v = av_y_v + y_v;
                            av_x_v = av_x_v + x_v;
                            
                            nr_vec = nr_vec + 1;
                    end  
                end
                av_vector_band_x{i}(j_b) = av_x_v / nr_vec;
                av_vector_band_y{i}(j_b) = av_y_v / nr_vec;
                av_activity_band{i}(j_b) = l_vec  / nr_vec;
            end
        end
        
     elseif operation == 7
        % correlate two vector maps by angle difference
        % read the vector map 1
        vecMap1 = load(fileList1{i}); 
        % read the vector map 2
        vecMap2 = load(fileList2{i});       
        
        
        a1=char(fieldnames(vecMap1));
        a2=char(fieldnames(vecMap2));
        for v=1:size(vecMap1.(a1),1)
            % scalar product
            v1_x = vecMap1.(a1)(v,3)- vecMap1.(a1)(v,1);
            v1_y = vecMap1.(a1)(v,4)- vecMap1.(a1)(v,2);
            v2_x = vecMap2.(a2)(v,3)- vecMap2.(a2)(v,1);
            v2_y = vecMap2.(a2)(v,4)- vecMap2.(a2)(v,2);
            l1 = sqrt(v1_x^2 + v1_y^2);
            l2 = sqrt(v2_x^2 + v2_y^2);
            scal_prod = (v1_x * v2_x +  v1_y * v2_y) / (l1 * l2);
            d_angle(nr,v,1) = vecMap1.(a1)(v,1);
            d_angle(nr,v,2) = vecMap1.(a1)(v,2);
            d_angle(nr,v,3) = scal_prod;
        end
        
        
    elseif operation == 8
        img1 = imread(fileList1{i}); 
        % corrected an normalized fret images as used for Rac
        res_image = ~img1;
    elseif operation == 9
        % subtract
        if vector_mode == 0
            img1 = imread(fileList1{i});
            if one_image_mode
                img2 = imread(fileList2{1});                
            else
                img2 = imread(fileList2{i});
            end
            res_image = double(img1) - double(img2);
        elseif vector_mode == 1
            img1 = imread(fileList1{i});
            res_image = double(img1) - av_intensity(i);
        end
        res_image = uint16((res_image > 0) .* res_image);
    elseif  operation == 10
        % subtract vector field
          % read the vector map 1
            vecMap1 = load(fileList1{i});
            % read the vector map 2
            vecMap2 = load(fileList2{i});

            a1=char(fieldnames(vecMap1));
            a2=char(fieldnames(vecMap2));
            for v=1:size(vecMap1.(a1),1)
                 dx = (vecMap1.(a1)(v,3)- vecMap1.(a1)(v,1)) - (vecMap2.(a2)(v,3)- vecMap2.(a2)(v,1));
                 dy = (vecMap1.(a1)(v,4)- vecMap1.(a1)(v,2)) - (vecMap2.(a2)(v,4)- vecMap2.(a2)(v,2));
                 dv(nr,v,1) = vecMap1.(a1)(v,1);
                 dv(nr,v,2) = vecMap1.(a1)(v,2);
                 dv(nr,v,3) = vecMap1.(a1)(v,1) + dx;
                 dv(nr,v,4) = vecMap1.(a1)(v,2) + dy;
                 %dl(round(vecMap1.(a1)(v,1)), round(vecMap1.(a1)(v,2))) = sqrt(dx^2+dy^2);
                 dl(nr,v,1) = vecMap1.(a1)(v,1); 
                 dl(nr,v,2) = vecMap1.(a1)(v,2);
                 dl(nr,v,3) = sqrt(dx^2+dy^2);
            end      
    elseif  operation == 11            
        % automatic image segmentation
        img1 = imread(fileList1{i});
        img_name = 'default';
        [ans, img_edge, res_image, edge_pixel, length_edge, frame_pos, rem_pix] = ...
                imFindCellEdge(img1, img_name, 0,'normalize',16);        
    elseif  operation == 13        
        % average intensity of non-zero pixels
        img1 = imread(fileList1{i});
        mask = imread(fileList2{i}); 
        img1_masked = double(img1) .* double(mask);
        nr_pix = sum(sum(mask));
        total_intensity = double(sum(sum(img1_masked)));
        av_intensity(i) = total_intensity / nr_pix;    
    elseif operation == 14    
        % get curvature from spline: 
        % k^2 = [dr^2 ddr^2 - (dr ddr)^2] / (dr^2)^3
        r_lower = x_spline(i).knots(1);
        r_upper = x_spline(i).knots(end); 
        p = r_lower:3:r_upper;
        x_spline_d  = fnder(x_spline(i));
        y_spline_d  = fnder(y_spline(i));
        x_spline_dd = fnder(x_spline(i),2);
        y_spline_dd = fnder(y_spline(i),2);
        
        %v_x    = fnval(x_spline(i), p);
        %v_y    = fnval(y_spline(i), p);
        v_x_d  = fnval(x_spline_d, p);
        v_y_d  = fnval(y_spline_d, p); 
        v_x_dd = fnval(x_spline_dd, p);
        v_y_dd = fnval(y_spline_dd, p); 
        
        t1 = (v_x_d  .* v_x_d    + v_y_d  .* v_y_d );
        t2 = (v_x_dd .* v_x_dd   + v_y_dd .* v_y_dd );        
        t3 = (v_x_d  .* v_x_dd   + v_y_d  .* v_y_dd) .^2; 
        t4 = t1.^3;
        kappa_square = (t1.*t2 - t3)./t4;
        kappa_square_int(i) = sum(kappa_square);
    elseif operation == 15
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%% Subtract background  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        img1 = imread(fileList1{i});
        if ~roi_bg            
            roi_mask = ~(imread(fileList2{i}));  
        end
        % get av bg intensity
        bg_img = double(roi_mask) .* double(img1);
        av_bg_int = sum(sum(bg_img))/sum(sum(roi_mask)); 
        % subtract
        res_image = uint16(double(img1)  -  av_bg_int);   
    elseif operation == 16 
        %translate
        img1 = imread(fileList1{i});
        xform = [ 1  0  0; 0  1  0; x_shift y_shift  1 ];
            
        tform_translate = maketform('affine',xform);
        res_image = imtransform(img1,tform_translate,...
                    'XData', [1 (size(img1,2))],...
                    'YData', [1 (size(img1,1))]); %,...
                    %'FillValues', .7,'bicubic'
        
        %res_image = uint16(double(img1));
    elseif operation == 17 
        % get translation parameters and correct image 2 
        img1 = imread(fileList1{i});
        img2 = imread(fileList2{i});
        
        
        if 1
            % resize image first to power of two
            [w h] = size(img1);
            crop_dim = min(w,h);
            crop_dim = 2^(nextpow2(max(size(img1)))-1);
            img1_crop = imcrop(img1,[h/2-crop_dim/2 w/2-crop_dim/2 crop_dim-1 crop_dim-1]);
            img2_crop = imcrop(img2,[h/2-crop_dim/2 w/2-crop_dim/2 crop_dim-1 crop_dim-1]);
    
            % filter image
            hsize = 20;
            sigma = 10;
            h = fspecial('gaussian',hsize,sigma);
            img1_crop_f = imfilter(img1_crop, h, 'replicate');
            img2_crop_f = imfilter(img2_crop, h, 'replicate');
            
            [phi(i), x_shift(i), y_shift(i)] = cytoAlignImage(img1_crop_f, img2_crop_f);
            res_image = imtransrot(img2, phi(i), x_shift(i), y_shift(i), 'bilinear');
        elseif 1 
             % filter image
            hsize = 20;
            sigma = 10;
            h = fspecial('gaussian',hsize,sigma);
            img1_f = imfilter(img1, h, 'replicate');
            img2_f = imfilter(img2, h, 'replicate');
            
            [res_image, x_shift(i), y_shift(i)] = kr_fretF(img1_f, img2_f);
            
            %res_image = imtransrot(img2, Phi, Row, Col, 'bilinear');
        else
            cytoGetImgShift(img1, img2);
        end
    elseif operation == 18
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%   Make movie  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % get colormap
        if strcmp(options.colormap,'hot')
            cmap = colormap(hot(256));
        elseif strcmp(options.colormap,'gray')
            cmap = colormap(gray(256));
        elseif strcmp(options.colormap,'jet')    
            cmap = colormap(jet(256));
        elseif strcmp(options.colormap,'hsv')
            cmap = colormap(hsv(256));
        elseif strcmp(options.colormap,'cool')
            cmap = colormap(cool(256));
        elseif strcmp(options.colormap,'summer')
            cmap = colormap(summer(256));    
        else
            cmap = color_map(grey(256));
        end
        
        img = imread(fileList1{i});

         
        if options.automatic_adjust 
            if ndims(img) < 3
                img = imadjust(img);
            end
            % convert image to double
            img = double(img);
            % down scale image to 8bit
            img = img ./ 2^ options.bit_depth .* 2^8;
        else
            % convert image to double
            img = double(img);
            
            % subtract offset
            img = img + options.offset;

            % down scale image to 8bit
            img = img ./ 2^ options.bit_depth .* 2^8;

            % stretch
            img = img .* options.scaling;
        end
        
        % convert image back to 8 bit
        img = uint8(img);
                        
        if options.img_resize ~= 1
            img = imresize(img,options.img_resize,'bicubic');
        end
              
        
        if options.bar
            % draw bar
            if strcmp(options.box_color,'black')
                box_color = 0;
            elseif  strcmp(options.box_color,'white')
                box_color = 255;
            end

            if ndims(img) < 3
                img(bar_coord(2)-2, bar_coord(1)-int8(bar_length/2):bar_coord(1)+int8(bar_length/2)) =  box_color;
                img(bar_coord(2)-1, bar_coord(1)-int8(bar_length/2):bar_coord(1)+int8(bar_length/2)) =  box_color;
                img(bar_coord(2),   bar_coord(1)-int8(bar_length/2):bar_coord(1)+int8(bar_length/2)) =  box_color;
                img(bar_coord(2)+1, bar_coord(1)-int8(bar_length/2):bar_coord(1)+int8(bar_length/2)) =  box_color;
            else
                img(bar_coord(2)-2, bar_coord(1)-int8(bar_length/2):bar_coord(1)+int8(bar_length/2),1:3) =  box_color;
                img(bar_coord(2)-1, bar_coord(1)-int8(bar_length/2):bar_coord(1)+int8(bar_length/2),1:3) =  box_color;
                img(bar_coord(2),   bar_coord(1)-int8(bar_length/2):bar_coord(1)+int8(bar_length/2),1:3) =  box_color;
                img(bar_coord(2)+1, bar_coord(1)-int8(bar_length/2):bar_coord(1)+int8(bar_length/2),1:3) =  box_color;
              
            end
        end

        if options.box 
            % draw box
            if strcmp(options.box_color,'black')
                box_color = 0;
            elseif  strcmp(options.box_color,'white')
                box_color = 255;
            end
            
            img(box_coord(1,2):box_coord(2,2), box_coord(1,1)) =  box_color;
            img(box_coord(1,2):box_coord(2,2), box_coord(2,1)) =  box_color;
            
            img(box_coord(1,2),  box_coord(1,1):box_coord(2,1)) =  box_color;
            img(box_coord(2,2),  box_coord(1,1):box_coord(2,1)) =  box_color;
        end
 
                   
        if options.timer 
            % add time stamp
            if ndims(img) < 3
                img(ind_text{i}) = uint8(255);
            else
                img(ind_text{i},1) = uint8(255);  
                img(ind_text{i},1) = uint8(255); 
                img(ind_text{i},1) = uint8(255); 
            end
        end
        
        f = im2frame(img,cmap);
        %text(20,20, time_string,'FontSize',30,'Color','w');
        avimovie = addframe(avimovie,f);
    elseif operation == 19
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%   Average images  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
        
         % average intensity of non-zero pixels
        img1 = imread(fileList1{i}); 
        mask = img1 > 0;
        counter_img = counter_img + mask;
        %nr_pix = sum(sum(mask));
        sum_img = sum_img+ double(img1);
    elseif operation == 20
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%   Overlay pixels  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
        img1 = imread(fileList1{i}); 
        % test if input image has double or 8 bit format
        if max(img1(:)) > 1 && max(img1(:)) <= 255
            edge_color = 255;
        elseif max(img1(:)) > 255
            edge_color = 2^16;
        else
            edge_color = 1;
        end
        % rainbow timer color coding
        cmap = colormap(jet(200)).*2^16;
        
        %img1 = 40.*cat(3,img1,img1,img1);
        
        if ndims(img1) == 2
            % BW images 
            %for i_img = 1 : i
            for i_pix = 1 : size(pixels{i},1)
                img1(pixels{i}(i_pix,2),pixels{i}(i_pix,1)) = edge_color;
                %img1(pixels{i_img}(i_pix,2),pixels{i_img}(i_pix,1)) = cmap(i,:);
                
                res_image = img1;
            end
            %end
        elseif ndims(img1) == 3
            % (rgb) color images
             for i_img = 1 : i
             for i_pix = 1 : size(pixels{i_img},1)
                %img1(pixels{i}(i_pix,2),pixels{i}(i_pix,1),1:3) = edge_color;
                img1(pixels{i_img}(i_pix,2),pixels{i_img}(i_pix,1),1:3) = cmap(i_img,:);
                res_image = img1;
             end
              end
        end
    elseif  operation == 21
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%   Average vector field   %%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
         vecMap1 = load(fileList1{i});
         
         a1=char(fieldnames(vecMap1));
         dx(i)=0;
         dy(i)=0;   
         vec_l(i) = 0;
         nr_vec = 0;

         in_poly = inpolygon(vecMap1.(a1)(:,2),vecMap1.(a1)(:,1), x_roi, y_roi);
         in_poly_vec = vecMap1.(a1)(find(in_poly),:);
         
         dx = in_poly_vec(:,3) - in_poly_vec(:,1);
         dy = in_poly_vec(:,4) - in_poly_vec(:,2);
         
         vec_l(i) = mean(sqrt(dx.^2+dy.^2));
         
         disp_l(i)= sqrt(mean(dx)^2+mean(dy)^2);
         
         % get histogram
         [vel_hist(i,:), vel_bins] = hist(sqrt(dx.^2+dy.^2),[0:0.3:7]);
    elseif  operation == 22
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%   Vectors from MPM matrix    %%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
         mpm_variable = load(fileList1{i});
         a=char(fieldnames(mpm_variable));
         if size(a,1) > 1
             MPM = mpm_variable.(a(2,:));
         else
             MPM = mpm_variable.(a);
         end
         mpm_to_vec(MPM, resultDir);
    elseif  operation == 23
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%  Analyse speckle track life time and speeds %%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
         mpm_variable = load(fileList1{i});
         a=char(fieldnames(mpm_variable));
         if size(a,1) > 1
             MPM = mpm_variable.(a(2,:));
         else
             MPM = mpm_variable.(a);
         end
         track_length_plotter(MPM, resultDir);  
    elseif  operation == 24
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%   Vectors from SCORES matrix    %%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
         scores_variable = load(fileList1{i});
         a=char(fieldnames(scores_variable));
         if size(a,1) > 1
             SCORES = scores_variable.(a(2,:));
         else
             SCORES = scores_variable.(a);
         end
         scores_to_vec(SCORES, resultDir);   
 	elseif  operation == 25
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%   Analyse protrusion         %%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
         protrusion_matrix = load(fileList1{i});
         a=char(fieldnames(protrusion_matrix));
         % (x,t)
         protrusion = protrusion_matrix.(a);     
         
         % segment protrusion
         p_m_p = protrusion > 0;
         Lp = bwlabel(p_m_p,4);
         events_p = regionprops(Lp,'Area');
         areas_p = [events_p.Area];
         % filter areas
         areas_p = areas_p(areas_p > 10);
         
         p_m_n = protrusion < 0;
         Ln = bwlabel(p_m_n,4);
         events_n = regionprops(Ln,'Area');
         areas_n = [events_n.Area];
         % filter areas
         areas_n = areas_n(areas_n > 10);
         
         av_areas_p = mean(areas_p);
         av_areas_n = mean(areas_n);
         area_ratio = av_areas_p / av_areas_n;
         
         % get average protrusion
         av_prot = mean(protrusion,2);
         
         
         nL = ['\n'];
         % get average protrusion
         av_av_prot = mean(av_prot);
         fprintf(1,nL); 
         disp('________________________________________________________');
         disp(fileList1{i});
         disp(['Average protrusion rate =                ', num2str(av_av_prot)]);
         
       
         % get std
         std_prot = std(protrusion,1,2);
         av_std = mean(std_prot);
         disp(['Average std of the protrusion rate =    ' , num2str(av_std)]);
         
         % get protrusion distribution
         protrusion_p = protrusion(protrusion > 0);
         protrusion_r = protrusion(protrusion < 0);
         
         protrusion_vec = reshape(protrusion,1,size(protrusion,1)*size(protrusion,2));
         p_vec = protrusion_vec(protrusion_vec  > 0);
         r_vec = protrusion_vec(protrusion_vec  < 0); 
         
         %[h,x]=hist(protrusion_vec,60);
         [h,x]=hist(protrusion_vec,60);
         % normalize
         %h = h./sum(h);
         h_area =  trapz(x,h);
         h = h./ h_area; 
         pos_indx = x > 0;
         neg_indx = x < 0;
         
         % protrusivity: ratio based on normalized histogram
         h_area_pos = trapz(x(pos_indx),h(pos_indx));
         h_area_neg = trapz(x(neg_indx),h(neg_indx)); 
         protrusivity_by_area_ratio = h_area_pos/h_area_neg;
         
         % protrusivity: difference based on normalized histogram
         protrusivity_by_area_diff = h_area_pos - h_area_neg;
         
         % edge movement rates
         av_prot_ret_rates = (mean(p_vec) - mean(r_vec))/2;
         
         
%          % protrusivity: ratio based on normalized histogram
%          net_protrusion = sum(h(pos_indx).*x(pos_indx));
%          net_retraction = sum(h(neg_indx).*x(neg_indx));
%          protrusivity_by_net = -net_protrusion/net_retraction;
         
         
         figure        
         plot(x,h);    
         title('Normalized protrusion, retraction velocity histogram');
         
         
         
         disp(['protrusivity by area ratio  =     ',  num2str(protrusivity_by_area_ratio)]);
         disp(['protrusivity by area diff   =     ',  num2str(protrusivity_by_area_diff)]);
         
         fprintf(1,nL); 
             
         disp(['Av protrusion area (pixel) =    ' , num2str(av_areas_p)]);
         disp(['Av retraction area (pixel) =    ' , num2str(av_areas_n)]);
         disp(['Ratio =    ' , num2str(area_ratio)]);
         
         fprintf(1,nL);
         
         disp(['Av prot ret events rates (um/min) =    ' , num2str(av_prot_ret_rates)]);
         
         
         %disp(['protrusivity by net advancement   =     ',  num2str(protrusivity_by_net)]);
         
         fprintf(1,nL);
         disp('                ');
         
         
         %protrusion_skewness = -skewness(protrusion_vec)
         %protrusion_kurtosis = kurtosis(protrusion_vec) 
         %protrusivity = protrusion_skewness * protrusion_kurtosis
%          x_model = 0:0.1:30;
%          model_curve = sin(x_model * 1);
%          figure
%          plot(x_model, model_curve);
%          [h_m,x_m]=hist(model_curve,60);
%          figure
%          plot(x_m, h_m);
         %prot_skew = skewness(protrusion_vec);   
         %prot_mean = mean(protrusion_vec);   
         %[h_prot,x_prot]=hist(protrusion_p,70);
         %[h_ret, x_ret]=hist(protrusion_r,70);
         %plot(x_prot, h_prot);
         %hold on
         %plot(x_ret, h_ret);

 	elseif  operation == 26
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%   Do PKA Fret                %%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
        if i == 1
            mkdir([resultDir filesep 'CFP_img_raw']);
            mkdir([resultDir filesep 'YFP_img_raw']);
        end
        
        img1 = imread(fileList1{i}); 
        %img2 = imread(fileList2);
        
        % left half is CFP; right half is YFP
        % get image size
        [img w img_h] = size(img1);
        
         % save split image
        num_str = sprintf('%4.4u', i)
        imwrite(uint16(img1(:,1:end/2)),     [resultDir filesep 'CFP_img_raw' filesep 'CFP_img_raw_' num_str '.tif'],'tif','Compression','none');
        imwrite(uint16(img1(:,end/2+1:end)), [resultDir filesep 'YFP_img_raw' filesep 'YFP_img_raw_' num_str '.tif'],'tif','Compression','none');
       
        % divide by flatfield image
        img = double(img1) ./ double(flatfield_img);
        
        % split image
        CFP_img = img(:,1:end/2);
        YFP_img = img(:,end/2+1:end);
        
        % corrected an normalized fret images as used for Rac
        if background 
            img1_c = double(img1) - av_bg1;
            img2_c = double(img2) - av_bg2;
        end
        
        % correct for bleedthrough
        % CFP to YFP bleedthrough
        %if YFP_correction
            YFP_img = YFP_img - 0.2 .* CFP_img;
        %end
        % YFP  to CFP bleedthrough
        % no correction 
        
        % calculation of the normalized and corrected FRET image
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        res_image = uint16((YFP_img ./ CFP_img) .* 5000);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % get rid of minus values outside of cell
        % res_image = res_image .* (res_image >= 0);
        % get rid of values > 1 outside of cell
        % res_image = res_image .* (res_image < 1);
        % res_image = uint16(res_image .* 2^16);
        % res_image = uint16(res_image .* 1000);
    end
    
    
    
    if operation == 3 || operation == 5 ||...
            operation == 8 || operation == 9 || operation == 15 ||...
            operation == 16 || operation == 17 || operation == 20 || operation == 26 
        [d1, currentFileName, currentFileExt] = fileparts(fileList1{i});
        imwrite(res_image, [resultDir filesep file_app currentFileName '.tif'],'tif','Compression','none');
    elseif operation == 11
        [d1, currentFileName, currentFileExt] = fileparts(fileList1{i});
        imwrite(res_image, [resultDir filesep file_app currentFileName '.tif'],'tif','Compression','none');   
        %[d1, currentFileName, currentFileExt] = fileparts(fileList1{i});
        %imwrite(img_edge, [resultDir filesep 'img_edge' currentFileName '.tif'],'tif','Compression','none');           
    elseif operation == 10
        %[d1, currentFileName, currentFileExt] = fileparts(fileList1{i});
        %save([resultDir filesep file_app currentFileName '_vel.mat'],'dv'); 
        %save([resultDir filesep file_app currentFileName '_speed.mat'],'dl'); 
    end
    fprintf(1,backSpc);
end
intwarning('on')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if operation == 1
    warning off all
    if stretch
        % stretch the multiplied image
        
        % get the maximal value
        max_total_img_val = max(max_img_val);
        % get the streching factor
        stretch_f = (2^14-1) / max_total_img_val;

        h_waitbar = waitbar(0,'Processing');
        for i=1:nr
            waitbar(i/nr, h_waitbar, [num2str(i)]);
            res_image_uint16 = uint16(stretch_f .* res_image{1,i});
            [d1, currentFileName, currentFileExt] = fileparts(fileList1{i});
            imwrite(res_image_uint16, [resultDir filesep file_app currentFileName '.tif'],'tif','Compression','none');
        end
    elseif stretch_i
        % stretch the multiplied image with individual factor
        
        % get the streching factor
        stretch_f = (2^14-1) ./ max_img_val;

        h_waitbar = waitbar(0,'Processing');
        for i=1:nr
            waitbar(i/nr, h_waitbar, [num2str(i)]);
            res_image_uint16 = uint16(stretch_f(i) .* res_image{1,i});
            [d1, currentFileName, currentFileExt] = fileparts(fileList1{i});
            imwrite(res_image_uint16, [resultDir filesep file_app currentFileName '.tif'],'tif','Compression','none');
        end        
    else
        h_waitbar = waitbar(0,'Processing');
        for i=1:nr
            waitbar(i/nr, h_waitbar, [num2str(i)]);
            res_image_uint16 = uint16(res_image{1,i});
            [d1, currentFileName, currentFileExt] = fileparts(fileList1{i});
            imwrite(res_image_uint16, [resultDir filesep file_app currentFileName '.tif'],'tif','Compression','none');
        end
    end
    clear res_image;
    close(h_waitbar);
    warning on    
end

if operation == 2
    % first show maximum value
    figure
    plot(max_img_val,'-x');
    xlabel('Frame');
    ylabel('Maximum intensity');
    % then write the maximum intensities into file
    dlmwrite([resultDir filesep file_app 'maximum_ratio_values.txt'], max_img_val);

    warning off all
    if stretch
        % stretch the divided image by one number
        % used that for time lapse series
        
        % get the maximal value
        max_total_img_val = max(max_img_val);
        % get the streching factor
        stretch_f = (2^14-1) / max_total_img_val;


        h_waitbar = waitbar(0,'Processing');
        for i=1:nr
            waitbar(i/nr, h_waitbar, [num2str(i)]);
            res_image_uint16 = uint16(stretch_f .* res_image{1,i});
            [d1, currentFileName, currentFileExt] = fileparts(fileList1{i});
            imwrite(res_image_uint16, [resultDir filesep file_app currentFileName '.tif'],'tif','Compression','none');
        end
    elseif stretch_i
        % stretch each individual divided image
        % get the maximal value

        % get the streching factor
        stretch_f = (2^14 - 1) ./ max_img_val;


        h_waitbar = waitbar(0,'Processing');
        for i=1:nr
            waitbar(i/nr, h_waitbar, [num2str(i)]);
            res_image_uint16 = uint16(stretch_f(i) .* res_image{1,i});
            [d1, currentFileName, currentFileExt] = fileparts(fileList1{i});
            imwrite(res_image_uint16, [resultDir filesep file_app currentFileName '.tif'],'tif','Compression','none');
        end       
    elseif stretch_c
        % stretch each image with user coefficient
        h_waitbar = waitbar(0,'Processing');
        for i=1:nr
            waitbar(i/nr, h_waitbar, [num2str(i)]);
            res_image_uint16 = uint16(stretch_coeff .* res_image{1,i});
            [d1, currentFileName, currentFileExt] = fileparts(fileList1{i});
            imwrite(res_image_uint16, [resultDir filesep file_app currentFileName '.tif'],'tif','Compression','none');
        end
    else
        h_waitbar = waitbar(0,'Processing');
        for i=1:nr
            waitbar(i/nr, h_waitbar, [num2str(i)]);
            res_image_uint16 = uint16(res_image{1,i});
            [d1, currentFileName, currentFileExt] = fileparts(fileList1{i});
            imwrite(res_image_uint16, [resultDir filesep file_app currentFileName '.tif'],'tif','Compression','none');
        end
    end
    clear res_image;
    close(h_waitbar);
    warning on
end

if operation == 4
    av_intensity_norm = av_intensity(1) ./ av_intensity;
    
    % do a quadratic fit of the bleach correction curve
    x = 1:length(av_intensity_norm);
    p = polyfit(x, av_intensity_norm, 2);
    av_intensity_norm_fit = p(1).* x.^2 + p(2) .* x + p(3);
    
    save([resultDir 'avImgIntensity.mat'],'av_intensity'); 
    save([resultDir 'integratedImgIntensityFit2.mat'],'av_intensity_norm_fit'); 
    
    figure
    plot(av_intensity_norm);
    hold on
    plot(av_intensity_norm_fit,'r');
    title('Correction factor (1/bleaching)');
end

if operation == 6
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%    ACTIVITY FROM EDGE      %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % convert cell to matrix
    if isempty(av_activity_band{1})
        av_activity_band(1) = [];
        nr = length(av_activity_band);
    end
    for i=1:nr
        l(i) = length(av_activity_band{i});
    end
    min_length = min(l);

    % create data matrix
    for i=1:nr
        av_activity_band_m(i,:) = av_activity_band{i}(1:min_length);
    end

    % transform pixel into metric length units
    if data.parameters.pixel_size_c
        x_axis = data.parameters.pixel_size.*(1:band_width:band_width*size(av_activity_band_m,2));
        xlabel_str = 'Distance from cell edge (microns)';
        if data.parameters.frame_interval_c
            ylabel_str = 'Velocity (microns/min)';
        else
            ylabel_str = 'Velocity (Pixel/frame)';            
        end
    else
        x_axis = 1:band_width:band_width*size(av_activity_band_m,2);
        xlabel_str = 'Distance from cell edge (pixel)';
        ylabel_str = 'Velocity (Pixel/frame)';
    end


    if scores_data
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%       Scores                %%%%%%%%%%%%%%%%%%%    
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
        if isempty(score_num{1})
            score_num(1) = [];
        end
        if isempty(score_density{1})
            score_density(1) = [];
        end        
        
        for i=1:nr
            score_num_m(i,:) = score_num{i}(1:min_length);
            score_density_m(i,:) = score_density{i}(1:min_length);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % get average scores per band
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %l = data.parameters.pixel_size;
%         if data.parameters.frame_interval_c
%             t = data.parameters.frame_interval.*60;
%             av_activity_band_m = av_activity_band_m ./ t;
%             ylabel_str = ('Average score value per band area and min');
%         else
            ylabel_str = ('Average score value per band area');
%       end
        
        av_scores = nanmean(av_activity_band_m,1);
        h=figure;plot(x_axis,av_scores);
        xlabel(xlabel_str);                
        ylabel(ylabel_str);
        print(h, [resultDir filesep 'av_scores_per_band.tif'],'-dtiff');
        hgsave(h, [resultDir filesep 'av_scores_per_band.fig']);        
        % save data
        av_scores_per_band = [x_axis', av_scores'];
        save([resultDir filesep 'av_scores_per_band'], 'av_scores_per_band');

        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % get the average number of scores per band
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        score_num_av = mean(score_num_m,1);
        h=figure;plot(x_axis,score_num_av);
        xlabel(xlabel_str);                
        ylabel('Number of scores per band');
        print(h, [resultDir filesep 'av_num_scores_per_band.tif'],'-dtiff');
        hgsave(h, [resultDir filesep 'av_num_scores_per_band.fig']); 
        % save data
        scores_number_per_band = [x_axis', score_num_av'];
        save([resultDir filesep 'scores_number_per_band'], 'scores_number_per_band');
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % get the average scores density per band
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if data.parameters.frame_interval_c
            score_density_m = score_density_m ./ t;
            ylabel_str = ('Scores density per band and min');
        else
            ylabel_str = ('Scores density per band');;
        end
        
        av_score_density = mean(score_density_m,1);
        h=figure;plot(x_axis,av_score_density);
        xlabel(xlabel_str);                
        ylabel(ylabel_str);
        print(h, [resultDir filesep 'scores_density_per_band.tif'],'-dtiff');
        hgsave(h, [resultDir filesep 'scores_density_per_band.fig']);         
        % save data
        density_per_band = [x_axis', av_score_density'];
        save([resultDir filesep 'density_per_band'], 'density_per_band');
        
        
    elseif images_data 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%       Images                %%%%%%%%%%%%%%%%%%%    
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
        % convert cell to matrix
        % create data matrix who's length is equal to the smalles cell
        % length
        for i=1:nr
            l(i) = length(av_activity_band{i});
        end
        min_length = min(l);

        for i=1:nr
            av_activity_band_m(i,:) = av_activity_band{i}(1:min_length);
        end       
        
        % plot the individual curves
        figure
        hold on
        for i=1:nr
            plot(av_activity_band_m(i,:));
        end
        
        % alignement of features 
        if 0
            
            
            
        end
        
        
        % ceate a image from the data matrix
        av_activity_band_m_resized = imresize(av_activity_band_m, 10, 'nearest');
        filter_kernel = fspecial('gaussian',[20 20], 3);
        av_activity_band_m_f = imfilter(av_activity_band_m_resized, filter_kernel, 'symmetric');
        figure,
        imshow(av_activity_band_m_f,[]);
        colormap(hot); colorbar
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % get intensity modulation in per band
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % get minimas and maximas
        if sum(sum((isnan(av_activity_band_m))))
            disp('No modulation analysis because of too few data points');   
        else
            h = fspecial('gaussian',7,4);
            for i=1:size(av_activity_band_m,2)
                av_activity_band_m_f2 = imfilter(av_activity_band_m(:,i),h,'replicate');
                max_p = find(imregionalmax(av_activity_band_m_f2));
                min_p = find(imregionalmin(av_activity_band_m_f2));
                % get ratio
                max_diff(i) = mean(av_activity_band_m_f2(max_p(1 : min(length(max_p),length(min_p))))-...
                    av_activity_band_m_f2(min_p(1 : min(length(max_p),length(min_p)))));
            end
            activity_modulation = max_diff./2;
            figure,plot(x_axis, activity_modulation);
            title('Modulation (average ampliture)')
            xlabel(xlabel_str);   
            ylabel('Modulation amplitude (intensity values)');
        end
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % get average intensity in band
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        av_av_activity_band_m   = nanmean(av_activity_band_m,1);
        std_av_activity_band_m  = nanstd(av_activity_band_m,1,1);
        if  normalise == 0
            av_av_n_activity_band_m = av_av_activity_band_m;
            std_av_n_activity_band_m  = std_av_activity_band_m;
        elseif  normalise == 1
            av_av_n_activity_band_m = av_av_activity_band_m ./ av_av_activity_band_m(1);
            std_av_n_activity_band_m  = std_av_activity_band_m ./ av_av_activity_band_m(1);
        elseif normalise == 2
            av_av_n_activity_band_m = av_av_activity_band_m ./ av_av_activity_band_m(end); 
            std_av_n_activity_band_m  = std_av_activity_band_m ./ av_av_activity_band_m(end);
        end
        
        
        figure,plot(x_axis, av_av_n_activity_band_m);
        hold on
        plot(x_axis, av_av_n_activity_band_m + std_av_n_activity_band_m,'--');
        plot(x_axis, av_av_n_activity_band_m - std_av_n_activity_band_m,'--');
        title('Mean and standart deviation');
        xlabel(xlabel_str);   
        if data.parameters.pixel_size_c
            excel_sheet_name = 'x_microns-y_activity';
        else       
            excel_sheet_name = 'x_pixel-y_activity';
        end
        ylabel('Activity');
        % save figure data as excel file
        xlswrite([resultDir 'data_from_edge'], [x_axis', av_av_n_activity_band_m'],excel_sheet_name);
        

        % save matrix
        %[pathstr, name, ext, versn] = filepartsfileparts(data.image1_path);
        save([resultDir 'av_activity_band_m.mat'],'av_activity_band_m');
        m_temp = [x_axis,activity_modulation];
        save([resultDir 'activity_modulation.mat'],'m_temp')
        %imwrite(res_image, [resultDir filesep file_app currentFileName '.tif'],'tif','Compression','none');
        
        % show band for visual control
        figure; imshow(res_image,[]);
    else 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%       Vectors               %%%%%%%%%%%%%%%%%%%    
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        for i=1:nr
            l(i) = length(av_activity_band{i});
        end
        min_length = min(l);

        %create data matrix
        for i=1:nr
            av_vector_band_x_m(i,:) = av_vector_band_x{i}(1:min_length);
            av_vector_band_y_m(i,:) = av_vector_band_y{i}(1:min_length);
            av_activity_band_m(i,:) = av_activity_band{i}(1:min_length);
        end
        % put it into one matrix so we can save it later
        av_vector_band_m = cat(3, av_vector_band_x_m, av_vector_band_y_m);
        
        % get speeds 
        av_speed_band = sqrt(av_vector_band_x_m.^2 + av_vector_band_y_m.^2);
        
        av_av_speed_band  = nanmean(av_speed_band,1);
        std_av_speed_band = nanstd(av_speed_band,1,1);
        
        if data.parameters.pixel_size_c && data.parameters.frame_interval_c
            av_av_speed_band = av_av_speed_band .* data.parameters.pixel_size ./...
            data.parameters.frame_interval.*60;
        end
        h=figure;plot(x_axis,  av_av_speed_band);
        hold on
        %plot(x_axis,  av_av_speed_band + std_av_speed_band, '--');
        %plot(x_axis,  av_av_speed_band - std_av_speed_band, '--');
        title('Average speed from average velocity');
        xlabel(xlabel_str);   
        ylabel(ylabel_str);
        print(h, [resultDir filesep 'av_speed_from_av_velocity.tif'],'-dtiff');
        hgsave(h, [resultDir filesep 'av_speed_from_av_velocity.fig']);  
        
        av_speed = nanmean(av_activity_band_m,1);
        if data.parameters.pixel_size_c && data.parameters.frame_interval_c
            av_speed = av_speed .* data.parameters.pixel_size ./...
            data.parameters.frame_interval.*60;
        end
        h=figure;plot(x_axis,  av_speed);
        title('Average speed');
        xlabel(xlabel_str);   
        ylabel(ylabel_str);  
        print(h, [resultDir filesep 'av_speed.tif'],'-dtiff');
        hgsave(h, [resultDir filesep 'av_speed.fig']);  
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % get modulation per band
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % check for Nan
        if sum(sum((isnan(av_speed_band))))
            disp('No modulation analysis because of too few data points');
        else
            % get minimas and maximas
            h = fspecial('gaussian',5,1);
            for i=1:size(av_speed_band,2)
                av_speed_band_f = imfilter(av_speed_band(:,i),h,'replicate');
                max_p = find(imregionalmax(av_speed_band_f));
                min_p = find(imregionalmin(av_speed_band_f));
                % get ratio
                max_diff(i) = mean(av_speed_band_f(max_p(1 : min(length(max_p),length(min_p))))-...
                    av_speed_band_f(min_p(1 : min(length(max_p),length(min_p)))));
            end
            figure,plot(x_axis, max_diff ./ 2);
            title('Modulation')
            xlabel(xlabel_str);   
            ylabel('Modulation amplitude (pixel)');
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % get average per band
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
%         av_speed_band_resized = imresize(av_speed_band, 10, 'nearest');
%         filter_kernel = fspecial('gaussian',[20 20], 3);
%         av_speed_band_resized_f = imfilter(av_speed_band_resized, filter_kernel, 'symmetric');
%         figure,imshow(av_speed_band_resized_f .* 6.4,[]);
%         axis on
%         % ticks
%         max_tick = size(av_vector_band_y_m, 2);
%         d_tick   = floor(max_tick / 10);
%         ticks = 0:dl:max_l;
%         if data.parameters.pixel_size_c
%             % labels 
%             for i = 1:length(ticks)
%                 max_l = size(av_vector_band_y_m, 2)*pixel_size;
%                 d_l   = floor(max_l / 10);
%             end
%         else
%             
%         end
%        set(gca,'XTicaxis onk',ticks)
%        set(gca,'XTickLabel',{'-pi','-pi/2','0','pi/2','pi'})
        colormap(jet);
        colorbar;
       
        
        figure,plot(x_axis,  av_av_speed_band);
        hold on
        plot(x_axis,  av_av_speed_band + std_av_speed_band, '--');
        plot(x_axis,  av_av_speed_band - std_av_speed_band, '--');
        title('Mean and standart deviation');
        if data.parameters.pixel_size_c
            xlabel('Distance from cell edge (micro m)');
        else       
            xlabel('Distance from cell edge (pixel)');
        end
        ylabel('Speed (pixel)');
        
         % save matrix
        %[pathstr, name, ext, versn] = fileparts(data.image1_path);
        %while pathstr(end-1),'\'
            
            
        save([resultDir 'av_speed_band_m.mat'],'av_speed_band');
        save([resultDir 'av_vector_band_m.mat'],'av_vector_band_m');
        set(handlesCytoProbe.info,'String',['Results saved to:  ', resultDir]);
        % show band for visual control
        figure; imshow(res_image,[]);

%         colormap(jet);
%         colorbar;
%        
%         
%         figure,plot(x_axis,  av_av_speed_band);
%         hold on
%         plot(x_axis,  av_av_speed_band + std_av_speed_band, '--');
%         plot(x_axis,  av_av_speed_band - std_av_speed_band, '--');
%         title('Mean and standart deviation');
%         xlabel(xlabel_str);   
%         ylabel('Speed (pixel)');
%         
%         % save matrix
%         save([resultDir 'av_speed_band_m.mat'],'av_speed_band');
%         save([resultDir 'av_vector_band_m.mat'],'av_vector_band_m');
%         set(handlesCytoProbe.info,'String',['Results saved to:  ', resultDir]);
%         % show band for visual control
%         figure; imshow(res_image,[]);
    end
end 

if operation == 8
    
    
end

if operation == 10
    % subtraction
    % average
    %if ndims(dl) > 2
    av_dl = squeeze(mean(dl,1));
    %end
    
    x=1:1:max(av_dl(:,1));
    y=1:1:max(av_dl(:,2));    
    [X,Y] = meshgrid(x,y);
    Z = griddata(av_dl(:,1),av_dl(:,2),av_dl(:,3),X,Y,'linear');    
    figure,surface(X,Y,Z);
    
    av_dv = squeeze(mean(dv,1));
    %VX = griddata(av_dv(:,1),av_dv(:,2),av_dv(:,3),X,Y,'linear');
    %VY = griddata(av_dv(:,1),av_dv(:,2),av_dv(:,4),X,Y,'linear');
    %figure,quiver(VX-X,VY-Y,10);
    figure,quiver(av_dv(:,1),av_dv(:,2),av_dv(:,3)-av_dv(:,1),av_dv(:,4)-av_dv(:,2),5);
end

if operation == 11
    % automatic segmentation
end

if operation == 13
    % bleach correction
    %av_intensity_norm = av_intensity(1) ./ av_intensity;
    
    % do a quadratic fit of the bleach correction curve
    l_x = length(av_intensity);
    x = 1:l_x;
    if 0
        % polynimial fit
        p = polyfit(x, av_intensity, 2);
        av_intensity_fit = p(1).* x.^2 + p(2) .* x + p(3);
    else
        % exponential fit
        [param resiudal] = fminsearch(@exp_fun, [av_intensity(1), log(av_intensity(end))/(l_x-1), 0.8],[], x, av_intensity);
        A = param(1);
        lambda = param(2);
        offset = param(3);
        av_intensity_fit = A .* exp(-lambda .* x) + offset;
    end
    
    av_intensity_fit_norm = av_intensity_fit(1) ./ av_intensity_fit;
     
    
    figure
    plot(av_intensity ./ av_intensity_fit(1));
    hold on
    plot(1 ./ av_intensity_fit_norm,'r');
    title('(Bleaching)');  
    
    % correct images
    h_waitbar = waitbar(0,'Processing');
    warning off all
    intwarning('off')
    for i=1:nr
        waitbar(i/nr, h_waitbar, [num2str(i)]);
        img1 = imread(fileList1{i});
        res_image = uint16(double(img1) * av_intensity_fit_norm(i));     
        % save images
        [d1, currentFileName, currentFileExt] = fileparts(fileList1{i});
        imwrite(res_image, [resultDir filesep file_app currentFileName '.tif'],'tif','Compression','none');    
    end
    set(handlesCytoProbe.info,'String',['Results saved to:  ', resultDir]);
    close(h_waitbar);
end

if operation == 14
    % curvature
    nr_el = size(kappa_square_int,2);
    mean_kappa_square_int = mean(kappa_square_int);
    std_kappa_square_int  = std(kappa_square_int);
    figure,plot(kappa_square_int);
    hold on
    plot(1:nr_el, mean_kappa_square_int .* ones(nr_el,1));
    plot(1:nr_el, (mean_kappa_square_int+std_kappa_square_int) .* ones(nr_el,1),'--');
    plot(1:nr_el, (mean_kappa_square_int-std_kappa_square_int) .* ones(nr_el,1),'--');
    title('Curvature square');
    xlabel('Frame #');
    ylabel('Curvature square');
end

if operation == 17
    % get transormation parameters
    av_x_shift  = mean(x_shift)
    av_y_shift  = mean(y_shift)
    
    if exist('phi')
        av_phi      = mean(phi)
        figure
        plot(phi);
        title('rotation (phi)');
    end
    figure
    plot(x_shift);
    hold on 
    plot(y_shift,'r')
    title('Translations');
end

if operation == 18
    % finish movie making
    avimovie = close(avimovie);
end

if operation == 19
    % calculate average image
    res_img1 =  sum_img ./ nr;
    res_img2 =  sum_img ./ counter_img;
    

    %imshow(res_img,[]);
    %imtool(res_img1);

    imtool(res_img2);
    % get histogram
    vec_img = reshape(res_img2,1,[]);
    vec_img_no = nonzeros(vec_img);
    [n,xout] = hist(vec_img_no,400);
    figure
    plot(xout,n);
end

if operation == 21
    % calculate average vector field

    figure,plot(vec_l);
    hold on
    plot(disp_l,'--');
    
    % get coherence measure 
    figure
    plot(disp_l ./ vec_l)
    title('Coherence');
    
    mean_hist = mean(vel_hist,1);
    figure,plot(vel_bins, mean_hist);
    title('Average speed histogram');
end


if operation == 1 || operation == 2 || operation == 3 || operation == 5 || operation == 6 ||...
   operation == 8 || operation == 9 || operation == 10|| operation == 11
    set(handlesCytoProbe.info,'String',['Results saved to:  ', resultDir]);
end

disp([num2str(nr) ' images processed']);
fclose all;
end

function residual = exp_fun(params, x, y)
        A = params(1);
        lambda = params(2);
        offset = params(3);
        FittedCurve = A .* exp(-lambda .* x) + offset;
        ErrorVector = FittedCurve - y;
        ErrorVector(1)=10*ErrorVector(1);
        residual = sum(ErrorVector .^ 2);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%














