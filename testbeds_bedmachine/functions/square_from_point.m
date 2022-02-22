function [X_crop, Y_crop, Q_crop] = square_from_point(X,Y, Q, data_i,n)
%SQUARE_FROM_POINT find the square of a given size with the datapoint at
%the center
%
%   Input:
%       X: X meshgrid of the original map
%       Y: Y meshgrid of the original map
%       datapoint: subscript indices of the given data point with in the original map
%       n: half number of points to sample at each direction
%          -> hence the side length is 2*n+1
%
%   Output:
%       
    x_i = data_i(2);
    y_i = data_i(1);
    leftmost = x_i - n;
    rightmost = x_i + n;
    upmost = y_i - n;
    downmost = y_i + n;
    
    % crop out by these vertices
    X_crop = X(upmost:downmost, leftmost:rightmost);
    Y_crop = Y(upmost:downmost, leftmost:rightmost);
    Q_crop = Q(upmost:downmost, leftmost:rightmost);
    
    
end

