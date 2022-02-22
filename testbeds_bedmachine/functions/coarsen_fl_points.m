function fl_data_coarse = coarsen_fl_points(fl_data, X, Y)
%COARSEN_FL_POINTS This function interpolate flow line data points at
%higher resolution to lower resolution by nearest point search to get
%non-repeated points in the lower resolution meshgrid, and we append to a
%new list
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

    % we are doing this one point at a time
    % if the dsearchn point has been found before, we do not keep it
    % otherwise, we append it to the list
    fl_data_coarse = [];
    for i = 1:size(fl_data, 1)
        temp_point_i = dsearchn([X(:), Y(:)], fl_data(i,:));
        % convert index to actual values
        temp_point   = [X(temp_point_i), Y(temp_point_i)];
        
        if isempty(fl_data_coarse)
            fl_data_coarse = temp_point;
        elseif ~ismember(temp_point, fl_data_coarse, 'rows') % not found before
            fl_data_coarse = [fl_data_coarse; temp_point];
        else
            continue % skip
        end
    end
end

