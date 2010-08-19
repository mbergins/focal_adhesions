function imOut = imtransrot(A, ang, row, col, method)

[so(1) so(2)] = size(A);
phi = ang * pi / 180; % Convert to radians

rotate = maketform('affine', [ cos(phi)  sin(phi)  0; ...
                              -sin(phi)  cos(phi)  0; ...
                              row        col      1 ]);
% Coordinates from center of A
hiA = (so - 1)/2;
loA = -hiA;

hiB = hiA;
loB = loA;
sn = so;

boxA = maketform('box', so, loA, hiA);
boxB = maketform('box', sn, loB, hiB);
T = maketform('composite', [fliptform(boxB), rotate, boxA]);

if strcmp(lower(method), 'bilinear')
    R = makeresampler('cubic', 'fill');
elseif strcmp(lower(method), 'bilinear')
    R = makeresampler('linear', 'fill');
else
    R = makeresampler('nearest', 'fill');
end
imOut = tformarray(A, T, R, [1 2], [1 2], sn, [], 0);
