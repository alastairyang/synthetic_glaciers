function thetas = point_derivative_angle(x, y, n)
%POINT_SLOPE Find the first derivative at every point along a line using
%a customized gradient scheme to counter the zigzag-y noise in data
%   
%   Input:
%       y: y-coordinate of the glacier flow line
%       n: number of points on each side to average. e.g., if n = 1, it is central
%          differnece (1 to the left and 1 to the right)
%
%   Output:
%       thetas: angle 
    
    % starting and ending points that need special treatement
    av_y_start = sum(y(2:n+1))/n;
    av_x_start = sum(x(2:n+1))/n;
    % forward diff
    grad_start = (av_y_start - y(1))/(av_x_start-x(1));
    
    % end
    av_y_end = sum(y(end-n:end-1))/n;
    av_x_end = sum(x(end-n:end-1))/n;
    % backward diff
    grad_end = (y(end) - av_y_end)/(x(end) - av_x_end);
    
    grad_start_n = repmat(grad_start, [n,1]);
    grad_end_n   = repmat(grad_end, [n,1]);
    
    % the first n points are now populated
    % The rest of y's
    % averaging the x and y coordinates of neighboring n points
    N = numel(x);
    grad = zeros(N, 1);
    for i = n+1:N-n-1
        left_av_y  = sum(y(i-n:i-1))/n;
        left_av_x  = sum(x(i-n:i-1))/n;
        right_av_y = sum(y(i+1:i+n))/n;
        right_av_x = sum(x(i+1:i+n))/n;
        % central difference
        % accounting for the scenario where the the line is straight up
        % hence delta_x = 0
        if right_av_x - left_av_x == 0
            grad(i) = NaN; % mark, substitute with 90 degree directly later
        end
        grad(i) = (right_av_y- left_av_y)/(right_av_x - left_av_x);
    end
    
    % add the start and end 
    grad(1:n) = grad_start_n;
    grad(end-n+1:end) = grad_end_n;
	
    % convert to angle
    % for vertical segments, we replace with 90 degree
    thetas = atan(grad);
    thetas(isnan(grad)) = pi/2;
end

