function cleaned_binary = remove_edge_adhesions(threshed_image,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Setup variables and parse command line
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

i_p = inputParser;
i_p.FunctionName = 'REMOVE_EDGE_ADHESIONS';

i_p.addRequired('threshed_image',@(x)isnumeric(x) || islogical(x));
i_p.addParamValue('binary_shift',0,@islogical);

i_p.parse(threshed_image,varargin{:});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Main Program
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

edge_border = ones(size(threshed_image));
if (not(any(strmatch('binary_shift',i_p.UsingDefaults))))
    edge_border = bwperim(i_p.Results.binary_shift);
else
    edge_border = bwperim(edge_border);
end

[row_indexes,col_indexes] = ind2sub(size(threshed_image), find(edge_border));
edge_adhesions = bwselect(threshed_image,col_indexes,row_indexes,4);

cleaned_binary = threshed_image & not(edge_adhesions);
