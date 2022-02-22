function [data_mean, data_std] = map_variance(data, time)
%MAP_VARIANCE Visualize the data variance among models at a given time
% 
%   example:  map_variance(dataset, 2015)
%
%   Input:
%       data: data from each group/model in a cell array
%       time: a number (double), indicating a year
%
%   Output:
%       no output. Returns a plot.

    % first, get the datetime
    assert(isa(time, 'double'), 'The time is not a double')
    dt_interest = datetime(time, 06, 30);
    % allow for ambiguity by including the previous year
    dt_interest_lastyr = datetime(time-1, 01, 01);
    
    % get the time and data in the dataset
    sz = numel(data);
    
    % initalize a datetime array
    darr = datetime(zeros(sz,1), 0, 0);
    
    for i = 1:sz
        % bracket the time
        bracketed_time_logical = data{i,1}.time > dt_interest_lastyr ...
                       & data{i,1}.time < dt_interest;
        % find the last one
        bracketed_time_idx = find(bracketed_time_logical == 1);
        last_bkt_idx = bracketed_time_idx(end);
        % add the time to the allocated datetime array
        darr(i) = data{i,1}.time(last_bkt_idx);
        
        % coordinates: compare the coordinates to previous ones
        if i == 1
            coor_X_old = data{1,1}.X;
            % initialize 3-D arrays for storing var data
            data_slices = zeros(size(coor_X_old, 1), size(coor_X_old, 2), sz);
        else
            coor_X_new = data{i,1}.X;
            % if the new coordinate is different from the last one, meaning
            % we have a different coord for this size
            % tolerance is 1 meter
            if ~isempty(find((coor_X_new - coor_X_old < 1) == 0, 1))
                error('Coordinates are different')
            end
        end
                
        % get the data at the given year
        data_slices(:,:,i) = squeeze(data{i,1}.var(:,:,last_bkt_idx));
    end
    
    % calculate the mean and variance along the 3rd dimension
    data_mean = mean(data_slices, 3);
    data_std = sqrt(var(data_slices, 0, 3));
    
    
    %% visualize
    % finding grid limit
    X = coor_X_old;
    % still take Y coordinates from the first dataset
    Y = data{1,1}.Y;
    X_max = max(X(:));
    X_min = min(X(:));
    Y_max = max(Y(:));
    Y_min = min(Y(:));
    
    ratio = (Y_max - Y_min)/(X_max - X_min);
    fig_width = 1200;
    figure('position',[100,100,fig_width, fig_width/2.5*ratio]);
    
    % first plot: mean
    subplot(1,2,1)
    s = pcolor(X, Y, data_mean);
    s.FaceColor = 'interp';
    s.EdgeColor = 'none';
    colorbar
    title('Mean')
    xlim([X_min X_max])
    ylim([Y_min Y_max])
    
    % second plot: variance
    subplot(1,2,2)
    s2 = pcolor(X, Y, data_std);
    s2.FaceColor = 'interp';
    s2.EdgeColor = 'none';
    colorbar
    title('Standard Deviation')
    xlim([X_min X_max])
    ylim([Y_min Y_max])
    

end

