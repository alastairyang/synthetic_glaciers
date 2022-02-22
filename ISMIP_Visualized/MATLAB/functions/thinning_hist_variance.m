function [] = thinning_hist_variance(fl_sector_data, model_data, data_for_plot)
%THINNING_HIST_VARIANCE Look for the agreement of timeseries' trends among
%models. This model will sample timeseries, find trends, test if it is
%ascending/descending trends, and aggregate into a value showing
%(dis)agreement.
%   
%   example:  thinning_hist_variance(flowline_coordinates, model_data)
%             (looking at only the SW sector of GrIS, all ISMIP model data)
%
%   Input:
%       fl_sector_data: a table(X,Y,Name) containing coordinates and labels
%       model_data: the cellarray containing ISMIP model data

    % Get timeseries cell array
    timeseries = point2timeseries(fl_sector_data, model_data, 1);
    
    n_model = numel(timeseries);
    n_fl    = numel(timeseries{1,1}.var);
    % create a cellarray storing all trend test stats
    trend_tests = cell(n_model,1);
    
    nan_count = 0;
        
    % iterate through models to 
    for i = 1:n_model
        % create Nx2 array storing decision for each data points
        trend_tests{i,1} = zeros(size(fl_sector_data, 1), 2);
        
        % iterate through flow lines
        for j = 1:n_fl
            ts_array = timeseries{i,1}.var{j,1}.data;
            n_ts_per_fl = size(ts_array, 2);
            % iterate through each data point on the flow line
            for k = 1:n_ts_per_fl
                % Hodrick-Prescott filter: decomposing into a trend and
                % cyclical component
                
                % first check if the timeseries is NaN
                if any(isnan(ts_array(:)))
                    nan_count = nan_count + 1;
                    sprintf('Found #%d matrix containing NaN', nan_count)
                    % convert all NaN to zeros
                    ts_array(isnan(ts_array) | isinf(ts_array)) = 0;
                end
               
                % get trend line
                trend = hpfilter(ts_array(:, k), 1600);
                % test if it is
                %   A: an ascending trend: F-stats p-value < .05, and
                %   coeff > 0
                %   B: a descending trend: F-stats p-value < .05, and coeff
                %   < 0
                % We are testing against an ascending trend
                
                X_asc = [transpose(1:numel(trend)), ones(numel(trend), 1)];
                [b, ~, ~, ~, stats] = regress(trend, X_asc);
                % look into stats_**, just look at the p_value (3rd one)
                if stats(3) <= 0.05 && b(1) > 0 
                    % ascend is significant
                    trend_class = 1;
                elseif stats(3) <= 0.05 && b(1) < 0
                    % descend is significant
                    trend_class = -1;
                else
                    % it is flat
                    trend_class = 0;
                end
                
                % assign the trend class to the output array cell
                index = (j-1)*n_ts_per_fl + k;
                trend_tests{i,1}(index, 1) = trend_class;
            end
        end
                
    end
    
    % now we have all the statistics 
    % combine all cells into the 2nd dimension of the array
    trend_tests_all = zeros(size(fl_sector_data, 1), n_model);
    for i = 1:n_model
        trend_tests_all(:,i) = trend_tests{i,1}(:,1);
    end
    
    % calculate the variance
    trend_std  = sqrt(var(trend_tests_all, 0, 2));
    trend_mean = mean(trend_tests_all, 2);
    
    % unpack fl_data into double arrays
    [coor, ~, ~] = unpack_fl_table(fl_sector_data);
    
    %% visualize
    % take the grid from the first model coordinates
    % finding grid limit
    X = model_data{1,1}.X;
    % still take Y coordinates from the first dataset
    Y = model_data{1,1}.Y;
    X_max = max(X(:));
    X_min = min(X(:));
    Y_max = max(Y(:));
    Y_min = min(Y(:));
    
    ratio = (Y_max - Y_min)/(X_max - X_min);
    fig_width = 1200;
    figure('position',[100,100,fig_width, fig_width/2.5*ratio]);
    
    % first plot: mean
    subplot(1,2,1)
    greenland();
    set(gca, 'visible','off')
    hold on
    %[zeros(numel(trend_std), 1), trend_mean, zeros(numel(trend_std), 1)]
    q1 = scatter(coor(:,1), coor(:,2), 15, trend_mean, 'filled');
%     q1.AlphaData = trend_mean;
%     q1.MarkerFaceAlpha = 'flat';
    colorbar
    title('Mean')
    xlim([X_min X_max])
    ylim([Y_min Y_max])
    axes('Position', [0.3, 0.18, 0.1, 0.15])
    set(gca, 'visible','off')
    box on
    histogram(trend_mean)
    hold off
    
    % second plot: variance
    subplot(1,2,2)
    greenland();
    set(gca, 'visible','off')
    hold on
    q2 = scatter(coor(:,1), coor(:,2), 15, trend_std, 'filled');
%     q2.AlphaData = trend_std;
%     q2.MarkerFaceAlpha = 'flat';
    colorbar
    title('Standard Deviation')
    xlim([X_min X_max])
    ylim([Y_min Y_max])
    axes('Position', [0.95-0.2, 0.18, 0.1, 0.15])
    box on
    histogram(trend_std)
    hold off
    
    %% histograms
    figure('position', [100, 100, 800, 300])
    subplot(1,2,1)
    histogram(trend_mean)
    title('Mean trend')
    
    subplot(1,2,2)
    histogram(trend_std)
    title('Trend standard deviation')
    
end

