function [] = visualize_flowline_bm(fl_data)
%VISULIZE_FLOWLINE_BM Visualizing flowlines which are for bedmachine
%extractions
%
%   Input:
%       fl_data[struc]: nested structures storing all flowline data
%
%   Output:
%       None
    
    glacier_names = fieldnames(fl_data);
    N_glacier = numel(glacier_names);
    for i = 1:N_glacier
        % number of flowlines
        flowline_names = fieldnames(fl_data.(glacier_names{i}));
        N_fl = numel(flowline_names);
        % choose the one close to/in the middle
        midnum_fl = middle_number(N_fl);
        % combine the glacier name and flowline name
        % and save to a new structure
        sampled_fl_data.([glacier_names{i}, '_', flowline_names{midnum_fl}])...
                  = fl_data.(glacier_names{i}).(flowline_names{midnum_fl});
    end
    
    % we obtain a structure of sampled flowline data
    
    %% Plot all data
    glacier_names = fieldnames(sampled_fl_data);
    N_glacier = numel(glacier_names);
    % get rid of "glacier" and "flowline..."
    % only keep the indices
    new_glacier_names = strings(N_glacier, 1);
    
    figure
    % plot greenland coastline
    greenland()
    
    for i = 1:N_glacier
        % first trim the names
        this_name  = glacier_names{i};
        start_trim = strfind(this_name, 'glacier');
        end_trim   = strfind(this_name, '_');
        new_glacier_names(i) = this_name(start_trim+numel('glacier'):end_trim-1);
        
        this_fl = sampled_fl_data.(glacier_names{i});
        hold on
        scatter(this_fl(:,1), this_fl(:,2), 'filled')
        % add the text to the first data point
        text(this_fl(1,1), this_fl(1,2) , new_glacier_names(i),...
            'VerticalAlignment','bottom','HorizontalAlignment','right')
    end
    hold off
    

end

