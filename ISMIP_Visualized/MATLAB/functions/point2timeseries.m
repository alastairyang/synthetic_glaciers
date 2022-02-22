function timeseries = point2timeseries(fl_data, dataset, normal)
%point2timeseries extract timeseries from given point coordinates
%
%   example: timeseries = point2timeseries(fl_data, dataset)
%
%   Input: 
%          fl_data: a table consisting of flow line data point coordiantes
%                   and their glacier name label
%          dataset: ISMIP datasets in a cell array
%          normal: 1 is normalization, 0 is no normalization
%
%   Output: 
%          timeseries cell array. Every cell is one ISMIP group/model
%          result...
%
%          In each cell, there are multiple cells representing
%          different glacier flow lines; within each of those cells we store
%          time series as vertical vectors, aggregating into a matrix
    
    % first determine if the X and Y of the datasets are 1-d array
    % or already a meshgrid
    N = numel(dataset);
    timeseries = cell(N,1);
    
    % unpacking fl_data from table to arrays
    [coor, names, unique_names] = unpack_fl_table(fl_data);
    n_fl = numel(unique_names);
    
    for i = 1:N
        data = dataset{i,1};
        if size(data.X) ~= size(data.Y)
            error('X and Y dimensions do not agree!')
        end
        
        % get meshgrid, if needed
        if any(size(data.X) == 1) % X, Y are in vector forms
            % meshgrid them
            [X, Y] = meshgrid(data.X, data.Y);
            if isinteger(X) || isinteger(Y)
                X = double(X);
                Y = double(Y);
            end
        else
            % if meshgrid by default, copy them
            X = data.X;
            Y = data.Y;
        end
        
        % get original meshgrid dimensions
        dimension = size(X);
        % reshape the meshgrid into two-column array for dsearchn 
        XY_vectorized = [reshape(X, [numel(X), 1]), ...
                            reshape(Y, [numel(Y), 1])];
        
        % look for the nearest points for one glacier flow line
        % Within each timeseries cell, create a cell array to store
        % timeseries at various points on the glacier
        timeseries_one_glacier = cell(n_fl, 1);
        time = data.time;
        n_t = numel(time);        
        % start to query each glacier
        for j = 1:n_fl
            % only get the corresponding glacier data points
            coor_select = coor(unique_names(j) == names, :);
            
            idx = dsearchn(XY_vectorized, coor_select);
            [row, col] = ind2sub(dimension, idx);
            % extract time series
            n_points = numel(row);
            
            % pre-allocate time matrix
            ts_all_points = zeros(n_t, n_points);
            
            % start to query each point on a glacier flow line
            for k = 1:n_points
                var_ts = data.var(row(k), col(k),:);
                % flatten the time object using reshape
                var_ts = reshape(var_ts, [n_t, 1]);
                % replace
                ts_all_points(:, k) = var_ts;
            end
            
            % NORMALIZATION
            if normal == 1
                % we consider the value at the first time step as reference
                % and diviving all values by that number as normalization
                ts_all_points = ts_all_points./ts_all_points(1,:);
            end
            
            % store all time series from that glacier
            %timeseries_one_glacier{j,1} = ts_all_points;
            % combine data and glacier name label into a struct
            each_glacier_fl.data = ts_all_points;
            each_glacier_fl.name = unique_names(j);
            
            % put this structure into a cell (it is not allowed to
            % concatenate multiple structs so this is necessary)
            timeseries_one_glacier{j} = each_glacier_fl;
        end
        
        % assign to a cell unit
        timeseries{i}.time = time;
        timeseries{i}.var = timeseries_one_glacier; % this is a cell
        timeseries{i}.name = data.att.institution;
        
        sprintf(' #%d file is completed', i)
    end
    
end

