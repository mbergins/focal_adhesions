








%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Low-pass image filtering
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The primary objective of this low pass filtering is really to prapare for
% segmentation and then, eventually, for centroid computation.
MedianWindSize = 9;   

for i = FrameMin:FrameMax 
    imgRaw = imread([handles.tempoutdir, 'RawTmp' num2str_fixwidth(i, 3) '.tif']);

    % Choose between median filter and Gauss filter
    imgMedian = medfilt2(imgRaw, [MedianWindSize MedianWindSize]);
    % imgMedian = Gauss2D(imgRaw, 7);
    
    imwrite(uint16(imgMedian), [handles.tempoutdir, 'medianTmp', num2str_fixwidth(i, 3), '.tif'], 'tif', 'Compression','none');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate Region
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
CropDim = 0;
StructElement = strel('disk', str2num(get(DilationTxt,'String')));
for i=FrameMin : FrameMax
    [LongDim, ThreshVal] = Threshold(tempoutdir, WindCalcs, i, get(handles.ThreshAbsRadio, 'Value'), str2num(get(handles.ThreshTxt,'String')),str2num(get(handles.ThreshSlopeTxt,'String')),StructElement);
    CropDim = max([CropDim LongDim]);
    arTransRot(i).ThreshVal = ThreshVal;
end

for i = FrameMin : FrameMax
    set(handles.FrameTxt,'String',['Resizing Frame ' num2str(i,'%1d')]);
    imgRaw    = imread([handles.tempoutdir, 'RawTmp' num2str_fixwidth(i, 3) '.tif']);
    imgMedian = imread([handles.tempoutdir, 'medianTmp' num2str_fixwidth(i, 3) '.tif']);
    rgnThresh = imread([handles.tempoutdir, 'rgnThreshTmp' num2str_fixwidth(i, 3) '.tif']);   % This temporary image has NOT be resized.
    rgnCentroid = regionprops(uint16(rgnThresh), 'Centroid'); % update to MATLAB v7
    
    % Center the images (based on computed centroid) 
    imgResize = imageResize(imgRaw, rgnCentroid.Centroid, handles.CropDim);
    imgMedianResize = imageResize(imgMedian, rgnCentroid.Centroid, handles.CropDim);
    rgnResize = imageResize(rgnThresh, rgnCentroid.Centroid, handles.CropDim);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % COMMENTS: Originally, the software uses the original image. This can
    % be a problem since the speckles are still there.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    imwrite(imgResize, [handles.tempoutdir, 'ResizeTmp' num2str_fixwidth(i, 3) '.tif'], 'tif', 'Compression', 'none');
    imwrite(imgMedianResize, [handles.tempoutdir, 'ResizeMedianTmp' num2str_fixwidth(i, 3) '.tif'], 'tif', 'Compression', 'none');
    imwrite(rgnResize, [handles.tempoutdir, 'rgnThreshResizeTmp' num2str_fixwidth(i, 3) '.tif'], 'tif', 'Compression', 'none');

    if(exist('thestack'))
        tempImg = imread(char(handles.chBImageStackList(i - handles.FrameMin + 1)));
        imgResize = imageResize(tempImg, rgnCentroid.Centroid, handles.CropDim);
        imwrite(imgResize,[handles.tempoutdir, '2ndWaveResizeTmp' num2str_fixwidth(i, 3) '.tif'],'tif','Compression','none');
    end
    
    % Basically, move the thresholded image such that its centroid is
    % aligned with the center of the expanded images 
    arTransRot(i).Trans1 = [handles.CropDim/2 handles.CropDim/2] - round(rgnCentroid.Centroid);  % computed translation
    arTransRot(i).Trans3 = [handles.CropDim/2 handles.CropDim/2] - round(rgnCentroid.Centroid);
    disp(['Resized Frame ' num2str(i) ', image size ' num2str(handles.CropDim) 'x' num2str(handles.CropDim) ' time ' num2str(toc)])
end




for image stack
    PhiRunning = 0;
    IRefFFT    = 0;

    [PhiRunning, Row, Col, IRefFFT] = AlignImage(IRefFFT, reference_image, translated_image, handles);
    
    
    RgnOut = imtransrot(RgnTst, PhiRunning, Row, Col, 'bilinear');
    ImgOut = imtransrot(ImgOri, PhiRunning, Row, Col, 'bilinear');
  

    arTransRot(i).Phi = PhiRunning;
    arTransRot(i).Trans2 = [Row Col];  
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function AlignImage
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [PhiRunningOut, RowOut, ColOut, IFFTOut] = AlignImage(IRefFFT, ImgTst, RgnTst, PhiRunningIn, handles)

%SEDisc=strel('disk',3);
%ImgTstOpen=imopen(ImgTst,SEDisc);
ImgTstOpen = ImgTst;
ImgTstMasked = double(ImgTstOpen)-mean(ImgTstOpen(logical(RgnTst))); % Set pixels outside the segmented region to be zero
ImgTstMasked(~logical(RgnTst))=0;  % set outside region to zero   

% So the function of RgnTst is to be a mask to set zero non-ROI pixels


if size(IRefFFT,1) == 1  % First time, when IRefFFT is ZERO.
    
    if strcmp(PhiRunningIn,'auto')
        RgnFeatures = regionprops(uint8(RgnTst),'Orientation');
        
        if (alignH ~= 0)
            PhiRunningIn = 0 - RgnFeatures.Orientation;  % Originally a 90 is added to make the spindle vertical. This is changed to 0 to make the spindle horizontal.
        elseif (alignV ~= 0)
            PhiRunningIn = 90 - RgnFeatures.Orientation;
        elseif (alignN ~= 0)
            PhiRunningIn = 0;
        end
    else
        PhiRunningIn = str2num(PhiRunningIn);  % set the angle offset directly  
    end
    
    PhiRunningOut = PhiRunningIn;
    RowOut = 0;
    ColOut = 0;
    
    IFFTOut = fft2(imrotate(ImgTstMasked, PhiRunningIn, 'bilinear', 'Crop'));
    % We need to compute this rotation to support offset
else
    [FitMaxRow, FitMaxCol, MaxInt] = CCMax(IRefFFT, ImgTstMasked, PhiRunningIn + 180, handles.CropDim);
    % MaxIntStats stores [Phi Row Col MaxVal]
    MaxIntStats = [0, FitMaxRow, FitMaxCol, MaxInt];  % Save the first run of maximum correlation. But at this point, no rotation yet (the first element, phi, is zero)
    
    PhiOptimal = 0;      % Current optimal value of phi
    MaxIntMax = MaxInt;  % Current best of correlation value
    
    figure(handles.WindCalcs);
    plot(MaxIntStats(1, 1), MaxIntStats(1, 4),'rx');
    %hold on;
    
    %%%%%%%%%%%%%%%%%%%%%%% START global searching %%%%%%%%%%%%%%%%%%%%%
    i=1; % counter of the number of stats   
    for PhiOffset = [1 -1]  % can only be either one of two directions
        phi = PhiOffset;
        FlagMax = 1;
        searchLimit = str2num(get(handles.PhiRangeTxt,'String'));
        while(FlagMax & abs(phi) <= searchLimit)
            i = i + 1;
            set(handles.FrameTxt,'String',['phi=' num2str(phi,'%1.2f')]);
            
            [FitMaxRow, FitMaxCol, MaxInt] = CCMax(IRefFFT, ImgTstMasked, PhiRunningIn + phi + 180, handles.CropDim);
            MaxIntStats(i,:) = [phi, FitMaxRow, FitMaxCol, MaxInt];
            
            FlagMax=0;
            if MaxInt > MaxIntMax
                FlagMax = 1;
                MaxIntMax = MaxInt;
                PhiOptimal = phi;
            end
            phi = phi + PhiOffset;  % Keep on searching in this direction
            
            figure(handles.WindCalcs);
            plot(MaxIntStats(:,1),MaxIntStats(:,4),'rx');
            drawnow;
        end
    end
    phi=PhiOptimal;
    %%%%%%%%%%%%%%%%%%%%%%% END global searching %%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%% START refined local searching %%%%%%%%%%%%%%%%%%%%%
    PowerRange = nextpow2(str2num(get(handles.RangePowerTxt, 'String'))) - 1;  % Get the precision defined by users
    for PowerCount = -1 : -1 : PowerRange
        for phi = [PhiOptimal - 2^PowerCount PhiOptimal+2^PowerCount]
            i = i + 1;
            set(handles.FrameTxt, 'String', ['phi=' num2str_fixwidth(phi,'%1.2f')]);
            
            [FitMaxRow, FitMaxCol, MaxInt] = CCMax(IRefFFT, ImgTstMasked, PhiRunningIn + phi + 180, handles.CropDim);
            MaxIntStats(i,:) = [phi, FitMaxRow, FitMaxCol, MaxInt];
            
            % Here MaxInt is really the maximum correlation 
            % The primary reason that the searching can be done is because
            % using FFT to compute correlation is fast !!
            
            if MaxInt > MaxIntMax
                MaxIntMax = MaxInt;
                PhiOptimal = phi;
            end
            
            figure(handles.WindCalcs);
            plot(MaxIntStats(:,1), MaxIntStats(:,4),'rx');
            drawnow;
        end
    end
    %hold off;
    %%%%%%%%%%%%%%%%%%%%%%% END refined local searching %%%%%%%%%%%%%%%%%%%%%
    %Fit best quadratic to last 6 values
    MaxIntStatsSort = sortrows(MaxIntStats, 4); % Sorting based on computed correlation value
    FitPhiMax = FindQuadMax(MaxIntStatsSort(end-5:end, 1),MaxIntStatsSort(end-5:end, 4)); % Find the best angle   
    
    PhiFitCloseness = abs(MaxIntStatsSort(:,1)-FitPhiMax);
    BestPhiIndex = find(PhiFitCloseness==min(PhiFitCloseness)); % Find the closest fit
    
    PhiRunningOut = PhiRunningIn + MaxIntStatsSort(BestPhiIndex,1);
    RowOut = MaxIntStatsSort(BestPhiIndex,2);
    ColOut = MaxIntStatsSort(BestPhiIndex,3);
    IFFTOut = fft2(imtransrot(ImgTstMasked, PhiRunningOut, RowOut, ColOut, 'bilinear'));  
    
    % Prepare for the plotting. Not really useful for optimization
    QuadCoeff = polyfit(MaxIntStatsSort(end-5:end,1),MaxIntStatsSort(end-5:end,4),2);
    QuadX = min(MaxIntStatsSort(end-7:end,1)) : (max(MaxIntStatsSort(end-7:end,1))-min(MaxIntStatsSort(end-7:end,1)))/50 : max(MaxIntStatsSort(end-7:end,1));
    QuadY = QuadCoeff(1)*QuadX.^2+QuadCoeff(2)*QuadX+QuadCoeff(3);
    plot(QuadX,QuadY, 'g-', MaxIntStatsSort(:,1), MaxIntStatsSort(:,4), 'rx', MaxIntStatsSort(end-5:end,1),MaxIntStatsSort(end-5:end,4),'bx', MaxIntStatsSort(BestPhiIndex,1),MaxIntStatsSort(BestPhiIndex,4),'mx');
    drawnow;
end




