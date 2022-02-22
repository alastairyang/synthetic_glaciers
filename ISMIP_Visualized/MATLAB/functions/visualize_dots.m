function [] = visualize_dots(data, coor, index)
%VISUALIZE Visualize the entire ice sheet in animation with given
%point coordinates plotted
%   Input: 
%       data - ice sheet dataset in a cell array; each cell a structure
%       coor - a table consisting of X and Y where X is longitude and Y is
%              latitude
%       index - the dataset of interest. Index query results from a
%       particular group/model
%   Output：
%       No output object. This function returns a plot.

    %% ice sheet data sets
    ds = data{index};
    try
        X = ds.X;
        Y = ds.Y;
    catch
        error('X,Y arrays not found')
    end
    
    % test if X,Y are meshgrid
    % first sample a vector
    X_v = ds.X(:,1);
    Y_v = ds.Y(1,:);
    % create a meshgrid from X_v, Y_v
    [X_test, Y_test] = meshgrid(X_v, Y_v);
    X_test = X_test';
    Y_test = Y_test';
    % compare if the same
    diff_X = abs(X_test - X);
    diff_Y = abs(Y_test - Y);
    % accounting for round-off error
    if max(diff_X(:)) > 1e-4 || max(diff_Y(:)) > 1e-4
        error('X, Y are not meshgird; not invariant along one dimension!')
    end
    
    % finding grid limit
    X_max = max(X(:));
    X_min = min(X(:));
    Y_max = max(Y(:));
    Y_min = min(Y(:));
    
    %% data points of interest
    % convert lats and longs from coor to polar stereographic
    lons = coor.X;
    lats = coor.Y;
    [X_p, Y_p] = ll2psn(lats, lons);
    coor_psn = [X_p, Y_p];
    
    % iterate over all coors and check if they are within the limits
    % keep a boolean index on the original table
    % (lying outside)
    XY_boolean = [X_p >= X_min, X_p <= X_max, Y_p >= Y_min, Y_p <= Y_max];
    % rows without zeros; along second dimension
    keep_row_boolean = all(XY_boolean, 2);
    
    % truncate 
    coor_psn_trim = coor_psn(keep_row_boolean, :);
    
    %% animate our plot
    ratio = (Y_max - Y_min)/(X_max - X_min);
    fig_width = 1200;
    fig1 = figure('position',[100,100,fig_width, fig_width/2.5*ratio]);
    subplot(1,2,1)
    % first, plot the data at t=1
    var_base = squeeze(ds.var(:,:,1));
    s = pcolor(X, Y, var_base);
    s.FaceColor = 'interp';
    s.EdgeColor = 'none';
    colorbar
    xlim([X_min X_max])
    ylim([Y_min Y_max])
    hold on
    scatter(coor_psn_trim(:,1), coor_psn_trim(:,2), 10, 'filled','r');
    hold off
    
    % the other subplot: animated
    % first, find the max and min values for the variables for all times
    for i = 2:numel(ds.time)
        var_diff = squeeze(ds.var(:,:,i)) - var_base;
        temp_min = min(var_diff(:));
        temp_max = max(var_diff(:)); 
        if i == 2 % first iteration, simply keep the values
            var_min = temp_min;
            var_max = temp_max;
            continue
        end
        % if more extreme values obtained, update
        if temp_min < var_min
            var_min = temp_min;
        end
        if temp_max > var_max
            var_max = temp_max;
        end
        
    end
    
    % with coloarbar limits obtained, visualize the animation
    % total animation time in seconds
    total_t = 5;
    dt = total_t/(numel(ds.time)-1);
    for i = 2:numel(ds.time)
        % first extract the variabla values;
        % we visualize the difference from the base (t=1)
        subplot(1,2,2)
        var_diff = squeeze(ds.var(:,:,i)) - var_base;
        s = pcolor(X, Y, var_diff);
        s.FaceColor = 'interp';
        s.EdgeColor = 'none';
        colorbar
        caxis([var_min var_max])
        xlim([X_min X_max])
        ylim([Y_min Y_max])
        t_str = datestr(ds.time(i));
        legend(t_str)
        hold on
        % plot datapoints, specifying markersize
        scatter(coor_psn_trim(:,1), coor_psn_trim(:,2),...
                10, 'filled','r');
        hold off
        pause(dt)

    end
    
end

