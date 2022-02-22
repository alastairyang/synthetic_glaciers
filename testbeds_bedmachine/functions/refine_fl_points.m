function fl_data_refine = refine_fl_points(fl_data, X, Y)
%REFINE_FL_POINTS interpolate the coarse line data points into higher
%resolution grid
%
%   Input:
%       fl_data[array]: Nx2 array containing the x's and y's of flowline
%                       datapoints.
%       X[array]:       meshgrid of x coordinates at lower resolution
%       Y[array]:       meshgird of y coordinates at lower resolution
%
%   Output:
%       fl_data_coarse[array]: Nx2 array containing x's and y's of flowline
%                       datapoints that are from the lower resolution map

    % Our strategy is to split the fl_data into multiple line segments
    % defined by two neighboring data points. We then interpolate wrt
    % BedMachine grid to get the data points in between and get their
    % nearest points on the grid.
    
    % we first get all the data points in fl_data on the grid
    defined_pts_i = dsearchn([X(:), Y(:)], fl_data);
    defined_pts   = [X(defined_pts_i), Y(defined_pts_i)];
    % remove repeated points
    defined_pts   = unique(defined_pts, 'row','stable');
    x = X(1,:);
    y = Y(:,1);
    
    N_pts = size(defined_pts, 1);
    N_pairs = N_pts - 1;
    fl_data_refine = [];
    for i = 1:N_pairs
        pt_1 = defined_pts(i  , :);
        pt_2 = defined_pts(i+1, :);
        pts  = [pt_1; pt_2];
        % interpolate with either (arbitrary) BedMachine axis
        % here we choose x axis
        [row1, col1] = ind2sub(size(X), defined_pts_i(i));
        x1_i = col1;
        [row2, col2] = ind2sub(size(X), defined_pts_i(i+1));
        x2_i = col2;
        if x1_i == x2_i % two points have the same x value, we should then use y axis
            % Mark: mark that we will eventually swap x and y to avoid
            % dealing with north-south lines that broke our coordinate
            % transformation system
            y1_i = row1;
            y2_i = row2;
            interv = 1;
            if y1_i > y2_i
                interv = -1;
            end
            this_y = y(y1_i:interv:y2_i);
            x_interp = interp1(pts(:,2), pts(:,1), this_y);
            this_y_betw = this_y(2:end-1); % exclude start and end
            this_x_betw = x_interp(2:end-1);
            
        else
            interv = 1;
            if x1_i > x2_i
                interv = -1;
            end
            this_x = x(x1_i:interv:x2_i);
            y_interp = interp1(pts(:,1), pts(:,2), this_x);
            this_x_betw = this_x(2:end-1);
            this_y_betw = y_interp(2:end-1);
        end
        % search the nearest in the XY grid for the interpolated pts
        betw_pts_i  = dsearchn([X(:), Y(:)], [this_x_betw(:), this_y_betw(:)]);
        betw_pts    = [X(betw_pts_i), Y(betw_pts_i)];
        % append to the refined flow line list
        new_pts = [pt_1; betw_pts];
        fl_data_refine = [fl_data_refine; new_pts];
        
        if i == N_pairs % append the last one
            fl_data_refine = [fl_data_refine; pt_2];
        end
    end
    
end

