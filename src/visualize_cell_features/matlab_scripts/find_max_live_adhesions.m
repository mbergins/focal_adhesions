function max_live_adhesions = find_max_live_adhesions(all_seqs)
% FIND_MAX_LIVE_ADHESIONS    determines the maximum number of adhesions
%                            tracked in each of the time steps in a given
%                            tracking matrix
%
%   M = find_max_live_adhesions(T_S) determines maximum number of
%   adhesions, 'M', tracked in a single time step from the tracking matrix,
%   'T_S'


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Setup variables and parse command line
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

i_p = inputParser;
i_p.FunctionName = 'FIND_MAX_LIVE_ADHESIONS';

i_p.addRequired('all_seqs',@isnumeric);

i_p.parse(all_seqs);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Main Program
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

max_live_adhesions = 0;

for i = 1:size(all_seqs,2)
    this_ad_num = length(find(all_seqs(:,i) > 0));
    
    if (this_ad_num > max_live_adhesions)
        max_live_adhesions = this_ad_num;
    end
end

end