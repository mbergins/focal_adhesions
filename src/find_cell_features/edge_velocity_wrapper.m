function edge_velocity_wrapper(varargin)

i_p = inputParser;
i_p.FunctionName = 'EDGE_VELOCITY_WRAPPER';

i_p.addOptional('file',@(x)exist(x,'file'));
i_p.addOptional('results',@ischar);
i_p.addOptional('max_img',@(x)isnumeric(x) && x >= 1);
i_p.addOptional('exclude_imgs',@ischar);
i_p.addOptional('contr',0,@(x)isnumeric(x) && (x == 1 || x == 0));
i_p.addOptional('protrusion',1,@(x)isnumeric(x) && (x == 1 || x == 0));
i_p.addOptional('t_step',1,@(x)isnumeric(x) && x >= 1);
i_p.addOptional('nojava',1,@(x)isnumeric(x) && (x == 1 || x == 0));

i_p.parse(varargin{:});

assert(isempty(strmatch('file',i_p.UsingDefaults)),'Error: must specify a starting file');
assert(isempty(strmatch('results',i_p.UsingDefaults)),'Error: must specify a folder for the results');
assert(isempty(strmatch('max_img',i_p.UsingDefaults)),'Error: must specify the maximum image number');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Main Program
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

exclude_imgs = [];
if (~isempty(i_p.Results.exclude_imgs))
    exclude_imgs = str2double(regexp(i_p.Results.exclude_imgs, ',', 'split'));
end

addpath(genpath('edge_velocity'));
if (i_p.Results.nojava)
    imEdgeTracker_nojava('contr',i_p.Results.contr,'protrusion',i_p.Results.protrusion,'t_step',i_p.Results.t_step,'file',i_p.Results.file,'results',i_p.Results.results,'max_img',i_p.Results.max_img,'exclude_imgs',exclude_imgs)
else
    imEdgeTracker('contr',i_p.Results.contr,'protrusion',i_p.Results.protrusion,'t_step',i_p.Results.t_step,'file',i_p.Results.file,'results',i_p.Results.results,'max_img',i_p.Results.max_img)
end