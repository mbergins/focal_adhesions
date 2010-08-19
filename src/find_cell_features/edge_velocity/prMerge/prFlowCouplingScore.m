function varargout = prFlowCouplingScore(VpFile,VcFile,shift)
%prFlowCouplingScore: This function reads two segment averaged flow-field files and outputs the 
%                     coupling score between the two flow fields. It can then be used for the analysis
%                     of the correlation to the edge protrusion.
%
% SYNOPSIS:
%    prFlowCouplingScore(VpFile,VcFile,scoreFile);
%    score = prFlowCouplingScore(VpFile,VcFile,scoreFile);
%
% INPUT:
%    VpFile    : The file that contains the primary flow average on segements.
%    VcFile    : The file that contains the coupling flow average on segements.
%    %scoreFile : The file where the coupling score is saved.
%    shift     : 0: save as 'score', 1: save as 'shift_score'
%
% OUTPUT:
%    score : (optional) A 2D matrix that contain the coupling score on segments and in time.
%            Row index corresponds to segment ID and column index corresponds to time steps.
%
% Author: Lin Ji
% Date  : Jan. 13, 2006

if nargout > 1
   error('Too many output arguments.');
end

%Load the two flow field.
S = load(VpFile);
a = fieldnames(S);
a_name = char(a{1});

Vp = S.(a_name);

S = load(VcFile);
a= fieldnames(S);
a_name = char(a{1});
Vc = S.(a_name);

score = sum(Vp.*Vc,3)./sum(Vp.^2,3);


if shift == 0
    save('score','score');
else
    shift_score = score;
    save('shift_score','shift_score');
end


if nargout == 1
   varargout{1} = score;
end
