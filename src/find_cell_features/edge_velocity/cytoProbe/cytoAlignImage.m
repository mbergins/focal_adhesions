function [PhiRunningOut, RowOut, ColOut] = cytoAlignImage(ImgTst, RgnTst)
% Original code by Ge Yang


ImgTstOpen = ImgTst;
ImgTstMasked = double(ImgTstOpen)-mean(ImgTstOpen(logical(ImgTst)));
% Set pixels outside the segmented region to be zero
ImgTstMasked(~logical(RgnTst))=0;  

RgnTstOpen = RgnTst;
RgnTstMasked = double(RgnTstOpen)-mean(RgnTstOpen(logical(RgnTst)));
% Set pixels outside the segmented region to be zero
RgnTstMasked(~logical(RgnTst))=0;  


FFT_ref = fft2(ImgTstMasked);

CropDim = 2^nextpow2(max(size(ImgTst)));
PhiRunningIn = 0;


[FitMaxRow, FitMaxCol, MaxInt] = CCMax(FFT_ref, RgnTstMasked, PhiRunningIn + 180, CropDim);
% Save the first run of maximum correlation. But at this point, no rotation
MaxIntStats = [0, FitMaxRow, FitMaxCol, MaxInt];  

PhiOptimal = 0;      % Current optimal value of phi
MaxIntMax = MaxInt;  % Current best of correlation value
    
%h_control=figure;
%plot(MaxIntStats(1, 1), MaxIntStats(1, 4),'rx');

    
%%%%%%%%%%%%%%%%%%%%%%% START global searching %%%%%%%%%%%%%%%%%%%%%
i=1; % counter of the number of stats
for PhiOffset = [1 -1]  % can only be either one of two directions
    phi = PhiOffset;
    FlagMax = 1;
    searchLimit = 10;
    while(FlagMax & abs(phi) <= searchLimit)
        i = i + 1;
        %disp(phi);

        %[FitMaxRow, FitMaxCol, MaxInt] = CCMax(FFT_ref, ImgTstMasked, PhiRunningIn + phi + 180, CropDim);
        [FitMaxRow, FitMaxCol, MaxInt] = CCMax(FFT_ref, RgnTstMasked, PhiRunningIn + phi + 180, CropDim);
        MaxIntStats(i,:) = [phi, FitMaxRow, FitMaxCol, MaxInt];

        FlagMax=0;
        if MaxInt > MaxIntMax
            FlagMax = 1;
            MaxIntMax = MaxInt;
            PhiOptimal = phi;
        end
        phi = phi + PhiOffset;  % Keep on searching in this direction

        %figure(h_control);
        %hold on;
        %plot(MaxIntStats(:,1),MaxIntStats(:,4),'rx');
    end
end
phi=PhiOptimal;
%%%%%%%%%%%%%%%%%%%%%%% END global searching %%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%% START refined local searching %%%%%%%%%%%%%%%%%%%%%
PowerRange = nextpow2(0.01) - 1;  % Get the precision defined by users
for PowerCount = -1 : -1 : PowerRange
    for phi = [PhiOptimal - 2^PowerCount PhiOptimal+2^PowerCount]
        i = i + 1;
        % disp(phi);

        %[FitMaxRow, FitMaxCol, MaxInt] = CCMax(FFT_ref, ImgTstMasked, PhiRunningIn + phi + 180, CropDim);
        [FitMaxRow, FitMaxCol, MaxInt] = CCMax(FFT_ref, RgnTstMasked, PhiRunningIn + phi + 180, CropDim);
        MaxIntStats(i,:) = [phi, FitMaxRow, FitMaxCol, MaxInt];

        % Here MaxInt is really the maximum correlation
        % The primary reason that the searching can be done is because
        % using FFT to compute correlation is fast !!

        if MaxInt > MaxIntMax
            MaxIntMax = MaxInt;
            PhiOptimal = phi;
        end

        %figure(h_control);
        %hold on;
        %plot(MaxIntStats(:,1), MaxIntStats(:,4),'rx');
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


% IFFTOut = fft2(imtransrot(ImgTstMasked, PhiRunningOut, RowOut, ColOut,
% 'bilinear'));
% Prepare for the plotting. Not really useful for optimization
%QuadCoeff = polyfit(MaxIntStatsSort(end-5:end,1),MaxIntStatsSort(end-5:end,4),2);
%QuadX = min(MaxIntStatsSort(end-7:end,1)) : (max(MaxIntStatsSort(end-7:end,1))-min(MaxIntStatsSort(end-7:end,1)))/50 : max(MaxIntStatsSort(end-7:end,1));
%QuadY = QuadCoeff(1)*QuadX.^2+QuadCoeff(2)*QuadX+QuadCoeff(3);
%plot(QuadX,QuadY, 'g-', MaxIntStatsSort(:,1), MaxIntStatsSort(:,4), 'rx', MaxIntStatsSort(end-5:end,1),MaxIntStatsSort(end-5:end,4),'bx', MaxIntStatsSort(BestPhiIndex,1),MaxIntStatsSort(BestPhiIndex,4),'mx');
%drawnow;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function CCMaxx
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [FitMaxRow, FitMaxCol, MaxInt] = CCMax(IRefFFTIn, ImgTstIn, phiIn, CropDim)   
ImgConv = fftshift(real(ifft2(IRefFFTIn.* fft2(imtransrot(ImgTstIn, phiIn, 0, 0, 'bilinear')))));
MaxInt = max(ImgConv(:));

% This is essentially computing the best translation parameters by using
% convolution computation. (Notice that the input images has already been
% translated to the center based on their calculated centroid

[MaxRow, MaxCol] = find(ImgConv == MaxInt);  

FitMaxRow = FindQuadMax([MaxRow-1 : MaxRow+1], ImgConv(MaxRow-1 : MaxRow+1, MaxCol)');  % Subpixel interpolation
FitMaxRow = FitMaxRow - CropDim/2;

FitMaxCol = FindQuadMax([MaxCol-1 : MaxCol+1], ImgConv(MaxRow,MaxCol-1 : MaxCol+1));    % Subpixel interpolation
FitMaxCol = FitMaxCol - CropDim/2;





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function FindQuadMax
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% f(x) = a * x^2 + b * x + c
% f'(x) = 2 * a * x + b = 0
% x = -b / (2 * a);
function xMax = FindQuadMax(x, y)
QuadCoeff = polyfit(x, y, 2);

if ~QuadCoeff(1)<0
    error('Error: Problem in quadratic curve fitting');
end

xMax = -(QuadCoeff(2)/(2*QuadCoeff(1)));


