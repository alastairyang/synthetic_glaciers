function [Xrot, Yrot] = rotation_wrt_point(X, Y, point, theta)
%%  Rotation of cartesian coordinates with respect to a specified point
%
%   Input:
%       X: meshgrid of X
%       Y: meshgrid of Y
%       point: coordinate (x,y) of the point that the rotation is with respect to
%       theta: rotational angle (starting the first quadrant
%
%   Output:
%       Xrot: rotated meshgrid of X
%       Yrot: rotated meshgrid of Y

    % displaced by the coordiantes of the point
    Xq = X - point(1);
    Yq = Y - point(2);

    Rot = [cosd(theta) -sind(theta); sind(theta) cosd(theta)];
    temp=[Xq(:),Yq(:)]*Rot' ;
          sz=size(Xq);

    Xrot=reshape(temp(:,1),sz);
    Yrot=reshape(temp(:,2),sz);

    % restoring the displaced distance
    Xrot = Xrot + point(1);
    Yrot = Yrot + point(2);
    
