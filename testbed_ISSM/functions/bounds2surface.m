function [new_data, new_X, new_Y] = bounds2surface(data, x, y, min_x, max_x, min_y, max_y)
%BOUNDS2SURFACE  This function takes in limits of both coordinates and
%output the patch of data bounded by the limits
%
%   example: bounds2surface(data, [1:100], [1:1000],...
%                           10, 20, 25, 35)
%            (x is from 10 to 20, y is from 25 to 35)
%
%   Input:
%       data: data of interest. Should be a meshgrid
%       x: x-coordinate in vector form
%       y: y-coordinate in vector form
%       min_x: lower bound of x
%       max_x: upper bound of x
%       min_y: lower bound of y
%       max_y: upper bound of y
%
%   Output:
%       new_data: the values of data in the patch
%       X: x-coordinates in meshgrid form
%       Y: y-coordinates in meshgrid form

    % find index and new meshgrid
    new_x_i = find(x<max_x & x>=min_x);
    new_y_j = find(y<max_y & y>=min_y);
    new_x = x(new_x_i);
    new_y = y(new_y_j);
    [new_X, new_Y] = meshgrid(new_x, new_y);

    % With the indices found, make a meshgrid
    % and transform to linear indexing
    [new_x_I, new_y_J] = meshgrid(new_x_i, new_y_j);
    sz = size(S);
    new_x_I_v = reshape(new_x_I, [numel(new_x_I), 1]);
    new_y_J_v = reshape(new_y_J, [numel(new_y_J), 1]);
    new_xy_linear = sub2ind(sz, new_y_J_v, new_x_I_v);

    % extract data and reshape to a meshgrid
    new_data_v = data(new_xy_linear);
    new_data = reshape(new_data_v, size(new_x_I));
    
    % visualize the patch
    meshgrid(new_X, new_Y, new_data)

end