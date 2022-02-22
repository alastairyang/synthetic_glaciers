function [X_new, Y_new, b, h, s] = straighten_flowline(x, y, map_X, map_Y, b, h, s, mk, n, dist_limit)
%%STRAIGHTEN_FLOWLINE Transform
%
%   Input:
%       x:     x-coordinate for the flow line
%       y:     y-coordinate for the flow line
%       map_X: BedMachine X mesh
%       map_Y: BedMachine Y mesh
%       b:     BedMachine base data
%       h:     BedMachine ice thickness data
%       s:     BedMachine surface elevation data
%       mk:    BedMachine mask data
%       n:     number of grid point -> half the side length for square
%       dist_limit: max length of flow line for our consideration, in km
%
%   Output:
%       X_new[array]: new X meshgrid for the straightened glacier geometry
%       Y_new[array]: new Y meshgrid for the straightened glacier geometry
%       b_interp_keep[array]: straightened glacier base
%       h_interp_keep[array]: straightened glacier ice thickness
%       s[array]:     straightened surface elevation
%
    %% parameters
    n_diff = 5; % numbers of point on each side in the numerical derivative
    
    %%
    gridsize = size(map_X);
    if size(x, 1) == 1
        xy = [x', y'];
    else
        xy = [x, y];
    end
    N_flpoints = size(xy, 1);

    % find angles for each point
    % _ points on one side
    thetas = point_derivative_angle(x, y, n_diff);

    b_interp_keep = [];
    h_interp_keep = [];
    s_interp_keep = [];
    total_dist = 0;
    
    for i = 1:N_flpoints
        % we only take every other point to avoid sharp corner, which is
        % the artifact from nearest point search
        % --> if even number, skip
         if rem(i, 2) == 0 
             continue
         else
            data = xy(i,:);
            % get the indices of this point
            [~, data_i_ind] = ismember(data, [map_X(:) map_Y(:)], 'rows' );
            [row, col] = ind2sub(gridsize, data_i_ind);
            % extract the square
            [X_square, Y_square, b_square] = square_from_point(map_X, map_Y, b, [row, col], n);
            [~       , ~       , h_square] = square_from_point(map_X, map_Y, h, [row, col], n);
            [~       , ~       , s_square] = square_from_point(map_X, map_Y, s, [row, col], n);
            % rotate the verticle line
            [Xrot, Yrot] = rotation_wrt_point(X_square, Y_square, data, thetas(i));
            % interpolate 
            b_interp = interp2(X_square, Y_square, b_square, Xrot, Yrot);
            h_interp = interp2(X_square, Y_square, h_square, Xrot, Yrot);
            s_interp = interp2(X_square, Y_square, s_square, Xrot, Yrot);
            % only keep the column parallel to y-axis and on the datapoint
            % append
            b_interp_keep = [b_interp_keep, b_interp(:,n+1)];
            h_interp_keep = [h_interp_keep, h_interp(:,n+1)];
            s_interp_keep = [s_interp_keep, s_interp(:,n+1)];

            % calculate flow line length, abort if two line
            % record current data point for calculation
            if i == 1 % just initiate
                last_data = xy(i,:);
            else 
                dist = sqrt((xy(i,1)-last_data(1,1))^2 +...
                            (xy(i,2)-last_data(1,2))^2);
                total_dist = (total_dist + dist)/1000; % convert to km
                if total_dist > dist_limit
                    % stop 
                    break;
                end

            end
        
        end
    end

    % make a meshgrid
%    [X_new, Y_new] = meshgrid(1:size(b_interp_keep, 1), 1:size(b_interp_keep, 2));
    
    % restore the geometry dimension to original grid
    total_dist = total_dist*1000; % convert to meter
    data_out = re_resol_alongflow(b_interp_keep, h_interp_keep, s_interp_keep, total_dist);
    b = data_out.b;
    h = data_out.h;
    s = data_out.s;
    X_new = data_out.X;
    Y_new = data_out.Y;
end