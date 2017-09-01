function [Img_out varargout] = atrous(Img_in,varargin)
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The algorithm proposed by Jean-Christophe in(2002, Pattern Recognitin)
% Extraction of spots in biological images using multiscale products
% usage: Img_out = atrous(Img_in), Img_out is grayscal without background
%        Img_out = atrous(Img_in,bg_flag,bw_flag);
%        Img_out = atrous(Img_in,bg_flag,bw_flag,exportLevel);
%        [Img_BW Img_BgFree] = atrous(Img_in);
%author: bei liu, 11/5/2009
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
switch nargin
    case 1
        bg_flag = 0;bw_flag = 0;th = 3;
    case 3
        bg_flag = varargin{1};bw_flag = varargin{2};th = 3;exportLevel = 1;
    case 4
        bg_flag = varargin{1};bw_flag = varargin{2};th = varargin{3};exportLevel = 1;
    case 5
        bg_flag = varargin{1};bw_flag = varargin{2};th = varargin{3};exportLevel = varargin{4};
    otherwise
        error('Invalid input parameters value!!!');
end

%progressText(O,'Perform atrous wavelet transform');
if iscell(Img_in)
    [row col] = size(Img_in{1});
    N = length(Img_in{1}); 
    temp = cell2mat(Img_in{1});
    clear Img_in
    Img_in = zeros(row,col,N);
    for i=1:N
        Img_in(:,:,i) = temp(((i-1)*row+1):i*row,:);
    end
else
    [row col N] = size(Img_in);    
end
Img_out = uint16(zeros(size(Img_in)));
switch nargout
    case 2
        varargout{1} = uint16(zeros(size(Img_in)));
end
J = 3;
A = cell(1,J);
W = cell(1,J);
%wavelet kernel
h = cell(1,J);
h{1} = [1/16 1/4 3/8 1/4 1/16];
h{2} = [1/16 0 1/4 0 3/8 0 1/4 0 1/16];
h{3} = [1/16 0 0 0 1/4 0 0 0 3/8 0 0 0 1/4 0 0 0 1/16];
if N>10
    multiWaitbar('Performing wavelet decomposition ...', 0);
    for i = 1:N
    % atrous wavelet transform
        for j = 1:J
            extn = (length(h{j})-1)/2;
            if j==1
                Temp = padarray(Img_in(:,:,i),[extn extn],'symmetric');
    %             Temp1 = conv2(h{j},h{j}',Temp,'same');
                Temp1 = imfilter(double(Temp),h{j},'symmetric','same');
                Temp1 = imfilter(Temp1,h{j}','symmetric','same');
                A{1,j} = Temp1(extn+1:extn+row,extn+1:extn+col);
                W{1,j} =  double(Img_in(:,:,i))-A{1,1};
            else
                Temp = padarray(A{1,j-1},[extn extn],'symmetric');
    %             Temp1 = conv2(h{j},h{j}',Temp,'same');
                Temp1 = imfilter(double(Temp),h{j},'symmetric','same');
                Temp1 = imfilter(Temp1,h{j}','symmetric','same');
                A{1,j} = Temp1(extn+1:extn+row,extn+1:extn+col);
                W{1,j} = A{1,j-1}-A{1,j};
            end
        end

        % %% Spot detection
        % %estimation of /sigma
        sigma = zeros(1,J); 
        t = zeros(1,J);
        for j= 1:J
            sigma(j) = (median(median(abs(W{1,j}-median(median(W{1,j}))))))/0.67;
            t(j) = th*sigma(j);
            W{1,j} = W{1,j}.*(W{1,j}>t(j));
    %        W{1,j} = (W{1,j}.^2-3*sigma(j)^2).*((W{1,j}.^2-3*sigma(j)^2)>0)./W{1,j};
        end
        % A_out = W{1,1}+W{1,2}+W{1,3}+A{1,3};
        if bw_flag == 1
            switch exportLevel
                case 1
                    Img_out(:,:,i) = sqrt(W{1,3}.*W{1,2})>1;
                case 2
                    Img_out(:,:,i) = 0.5*(W{1,3}+W{1,2})>1;
                case 3
                    Img_out(:,:,i) = W{1,2}>1;
                case 4
                    Img_out(:,:,i) = W{1,3}>1;
            end
        else
    %         Img_out(:,:,i) = sqrt(W{1,3}.*W{1,2});
            Img_out(:,:,i) = W{1,1}+W{1,2};
            if bg_flag == 1
                Img_out(:,:,i) = Img_out(:,:,i)+W{1,3};
            end
        end
        multiWaitbar('Performing wavelet decomposition ...', i/N);
       % progressText(i/N,'Perform atrous wavelet transform');
    end
    multiWaitbar('Performing wavelet decomposition ...', 'close');
else
    for i = 1:N
    % atrous wavelet transform
        for j = 1:J
            extn = (length(h{j})-1)/2;
            if j==1
                Temp = padarray(Img_in(:,:,i),[extn extn],'symmetric');
    %             Temp1 = conv2(h{j},h{j}',Temp,'same');
                Temp1 = imfilter(double(Temp),h{j},'symmetric','same');
                Temp1 = imfilter(Temp1,h{j}','symmetric','same');
                A{1,j} = Temp1(extn+1:extn+row,extn+1:extn+col);
                W{1,j} =  double(Img_in(:,:,i))-A{1,1};
            else
                Temp = padarray(A{1,j-1},[extn extn],'symmetric');
    %             Temp1 = conv2(h{j},h{j}',Temp,'same');
                Temp1 = imfilter(double(Temp),h{j},'symmetric','same');
                Temp1 = imfilter(Temp1,h{j}','symmetric','same');
                A{1,j} = Temp1(extn+1:extn+row,extn+1:extn+col);
                W{1,j} = A{1,j-1}-A{1,j};
            end
        end

        % %% Spot detection
        % %estimation of /sigma
        sigma = zeros(1,J); 
        t = zeros(1,J);
        for j= 1:J
            sigma(j) = (median(median(abs(W{1,j}-median(median(W{1,j}))))))/0.67;
            t(j) = th*sigma(j);
            W{1,j} = W{1,j}.*(W{1,j}>t(j));
    %        W{1,j} = (W{1,j}.^2-3*sigma(j)^2).*((W{1,j}.^2-3*sigma(j)^2)>0)./W{1,j};
        end
        % A_out = W{1,1}+W{1,2}+W{1,3}+A{1,3};
        if bw_flag == 1
            switch exportLevel
                case 1
                    Img_out(:,:,i) = sqrt(W{1,3}.*W{1,2})>1;
                case 2
                    Img_out(:,:,i) = 0.5*(W{1,3}+W{1,2})>1;
                case 3
                    Img_out(:,:,i) = W{1,2}>1;
                case 4
                    Img_out(:,:,i) = W{1,3}>1;
            end
        else
    %         Img_out(:,:,i) = sqrt(W{1,3}.*W{1,2});
            Img_out(:,:,i) = W{1,1}+W{1,2};
            if bg_flag == 1
                Img_out(:,:,i) = double(Img_out(:,:,i))+W{1,3};
            end
        end
    end
end