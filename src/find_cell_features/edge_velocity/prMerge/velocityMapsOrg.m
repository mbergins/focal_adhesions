function vel_field  = velocityMapsOrg(imgSize,M, MPM, time, mask_img, d0)

% Extract vectors
Mp=M(:,:,time);
Mv=Mp(find(Mp(:,1)~=0 & Mp(:,3)~=0),:);

% Mask non cell values
i=1;
while i <= size(Mv,1)
    if mask_img(round(Mv(i,1)), round(Mv(i,2))) == 0
        Mv(i,:)=[];
        i=i-1;
    end
    i=i+1;
end

if 1
    % Interpolate
    Md=vectorFieldInterpN(Mv,d0);
    vel_field = [Mv(:,1), Mv(:,2), Md];
else
    % just get the velocities
    Md = [Mv(:,3) - Mv(:,1), Mv(:,4) - Mv(:,2)];
    vel_field = [Mv(:,1), Mv(:,2), Md];
end