function [x_matrix, y_matrix] = extract_flowline_patch(x, y, n, ds)
%EXTRACT_FLOWLINE_PATCH extract coordinates of neighboring points
%   
%   Input:
%       x: a vector of x-coordinates of the flow line
%       y: a vector of y-coordinates of the flow line
%       n: number of neighboring points to extract on one side
%       ds: grid interval between extracted neighboing
%       points
%
%   Output:
%       x_matrix: basically x mesh
%       y_matrix: y coordinates of each extracted point in the patch
%       
    N = numel(x); % same as number of y
    y_matrix = zeros(2*n+1, N);
    
    for i = 1:N
        
        y_up = flip(y(i):-ds:y(i)-ds*n);
        y_down = y(i):ds:y(i)+ds*n;
        % y(i) is repeated twice; truncate one
        y_up = y_up(1:end-1);
        y_matrix(:,i) = [y_up, y_down]';
        
    end
    
    % repeat x vector into a matrix
    if size(x, 2) == 1
        x = x';
    end
    x_matrix = repmat(x, [2*n+1, 1]);
    
    
end

