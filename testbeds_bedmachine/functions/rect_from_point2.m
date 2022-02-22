function [map_X_crop, map_Y_crop, map_b_crop, map_h_crop, map_s_crop, map_mk_crop, ...
          map_X_cc, map_Y_cc, map_b_cc, map_h_cc, map_s_cc, map_mk_cc] =...
          rect_from_point2(fl_xy, map_X, map_Y, map_b, map_h, map_s, map_mk, n, dist_limit)
%RECT_FROM_POINT2 extracting a larger rectangle enclosing the flow line so
%that we can work with this instead of the entire BedMachine data
%   
%   Input:
%       fl_xy_i: a Nx2 array of flow line coordinates 
%       map_X:   X meshgrid of BedMachine
%       map_Y:   Y meshgrid of BedMachine
%       map_b:   bed topography map from BedMachine
%       map_h:   ice thickness map from BedMachine
%       map_mk:  ice/ocean mask map from BedMachine
%       n:       number of points as extension (larger than the minimal
%                bounding rectangle)
%       dist:    flow line distance limit
%   
%   Output:
%       map_*_crop: cropped full flow line length (rectangle shape) data
%       map_*_cc:   data further cropped to the provided flow line length
%                   limit
%       
    % fisrt, look for the coordinate bounds of the flow line
    x_max = max(fl_xy(:,1));
    x_min = min(fl_xy(:,1));
    y_max = max(fl_xy(:,2));
    y_min = min(fl_xy(:,2));
    
    % get map_x, map_y from map_X and map_Y (vector from meshgrid matrix)
    map_x = transpose(map_X(1,:));
    map_y = map_Y(:,1);
    
    % get the subscript for each dimension
    % add the padding
    x_max_i = dsearchn(map_x, x_max)+n;
    x_min_i = dsearchn(map_x, x_min)-n;
    % since the coordinate in matrix space along y-dimension is different
    % from the cartesian coordinate, so we minus n at y_max_i
    y_max_i = dsearchn(map_y, y_max)-n;
    y_min_i = dsearchn(map_y, y_min)+n;
    
    % finally, crop out the rectangle we want
    map_X_crop = map_X(y_max_i:y_min_i, x_min_i:x_max_i);
    map_Y_crop = map_Y(y_max_i:y_min_i, x_min_i:x_max_i);
    map_b_crop = map_b(y_max_i:y_min_i, x_min_i:x_max_i);
    map_h_crop = map_h(y_max_i:y_min_i, x_min_i:x_max_i);
    map_s_crop = map_s(y_max_i:y_min_i, x_min_i:x_max_i);
    map_mk_crop = map_mk(y_max_i:y_min_i, x_min_i:x_max_i);
    
    %% crop a smaller area based on the limit of flowline distance
    [~, index] = dist_coor_along_fl(fl_xy, dist_limit);
    fl_limit = fl_xy(1:index,:);
    
    % fisrt, look for the coordinate bounds of the flow line
    x_max = max(fl_limit(:,1));
    x_min = min(fl_limit(:,1));
    y_max = max(fl_limit(:,2));
    y_min = min(fl_limit(:,2));
    
    % get the subscript for each dimension
    % add the padding, 20 is enough
    x_max_i = dsearchn(map_x, x_max)+20; 
    x_min_i = dsearchn(map_x, x_min)-20;
    % since the coordinate in matrix space along y-dimension is different
    % from the cartesian coordinate, so we minus n at y_max_i
    y_max_i = dsearchn(map_y, y_max)-20;
    y_min_i = dsearchn(map_y, y_min)+20;
    
    % finally, crop out the rectangle we want
    map_X_cc = map_X(y_max_i:y_min_i, x_min_i:x_max_i);
    map_Y_cc = map_Y(y_max_i:y_min_i, x_min_i:x_max_i);
    map_b_cc = map_b(y_max_i:y_min_i, x_min_i:x_max_i);
    map_h_cc = map_h(y_max_i:y_min_i, x_min_i:x_max_i);
    map_s_cc = map_s(y_max_i:y_min_i, x_min_i:x_max_i);
    map_mk_cc = map_mk(y_max_i:y_min_i, x_min_i:x_max_i);
    
end

