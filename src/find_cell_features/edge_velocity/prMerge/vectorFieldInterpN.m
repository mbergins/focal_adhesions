function Mi=vectorFieldInterp(M,d0)
% vectorFieldInterp interpolates a vector field on a user-specified grid
%
% For a vector field v:  
%    Interpolation with a correlation matrix K : <v> = K * v      (convolution)
%    where:  K = sG*exp(-(dx^2+dy^2)/d0^2), sG = weigth vector
%
% SYNOPSIS   Mi=vectorFieldInterp(M,Pg,d0)
%
% INPUT      M       : vector field, stored in a (nx4)-matrix of the form [y0 x0 y x]n
%                      (where (y0,x0) is the base and (y,x) is the tip of
%                      the vector).
%            Pg      : regular grid points, stored in a (mx2)-matrix of the form [yg xg]m.
%            d0      : parameter for the weight function G=exp(-D.^2/(1+d0^2)),
%                      where D is the distance matrix between all grid
%                      points and all vector (base) positions.
%                      d0 can be a scalar or a matrix with size (nxm). See
%                      calcD0fromDiv.
%
% OUTPUT     Mi      : interpolated vector field.
%
% DEPENDENCES          vectorFieldInterp uses { createDistanceMatrix (C-MEX function) }
%                      vectorFieldInterp is used by { }
%
% Aaron Ponti, 11/18/2002

% Vector base positions
Pi=M(:,1:2);

% Vectors
V=[M(:,3)-M(:,1) M(:,4)-M(:,2)];

% Calculate distances
D=createDistanceMatrix(Pi,Pi);

% Correlation matrix (d0 may be a scalar or a matrix)
G=exp(-D.^2./d0.^2); clear D;

% Interpolate
Vi=[G*V(:,1) G*V(:,2)];

% Normalize
sG=sum(G,2);
sG(find(sG==0))=1; % Prevent division by zero
Vi=[Vi(:,1)./sG Vi(:,2)./sG];

% Mi is the interpolated M
%Mi=[Pg Pg+Vi];
Mi=Vi;