function [glacier_idx_parsed, nc_data] = read_felikson_nc(indices, dir)
%READ_FELIKSON_NC This file reads all flow line coordinate data (x,y) in NetCDF from
%Felikson et al., 2020
%
%   Input:
%       indices: indices of the glacier flow line
%       dir: directory to the folder containing all nc files
%
%   Output:
%       idx_parsed: parsed indices by glacier; stored as cell arrays within
%       strutures fields
%       nc_data: containing all extracted data

%% First parse the indices
% Group indices by nc file, since each flow line file can contain multiple
% flowlines XXX_01, XXX_02...; there also exist multiple branches for certain glaciers
% are denoted by character (e.g., a,b,c...)

 % first find if the indices contain "_"
 % if not, it is the glacier index followed by individual flow line of the
 % glacier
 N_idx = numel(indices);
 for i = 1:N_idx
     % check if it is a glacier main index
     if ~ismember('_', indices{i})
        main_idx = indices{i};
        full_idx = ['glacier',main_idx];
        % current number of cells
        glacier_idx_parsed.(full_idx) = {};
        N_cell = 1;
        % then it goes to the next iteration
     else
        % then looking for all flowlines that belong to this glacier
        % which is the previous main_idx
        if all(main_idx == indices{i}(1:numel(main_idx)))
            flowline_number = indices{i}(numel(main_idx)+2:end);
            % create a string that is consistent with the name in NetCDF
            fl_number_name = ['flowline', flowline_number];
            % add to the structure
            glacier_idx_parsed.(full_idx){N_cell} = fl_number_name;
            N_cell = N_cell+1;
        end
     end
     
 end
 
%% now that everything is stored 
% we extract nc data
glacier_filenames = fieldnames(glacier_idx_parsed);
N_glaciers = numel(glacier_filenames);
for i = 1:N_glaciers
    if isa(glacier_filenames{i}, 'string')
        this_glacier = convertStringsToChars(glacier_filenames{i});
    else
        this_glacier = glacier_filenames{i};
    end
    full_dir = [dir, this_glacier,'.nc']; 
    % inner loop: for each flow line of the glacier
    % get the number of flow line first
    N_flowline = numel(glacier_idx_parsed.(this_glacier));
    for j = 1:N_flowline
        if isa(glacier_idx_parsed.(this_glacier){j}, 'string')
            this_flowline = convertStringsToChars(glacier_idx_parsed.(this_glacier){j});
        else
            this_flowline = glacier_idx_parsed.(this_glacier){j};
        end
        % first try, since it seems that not every flow line index has data
        try
            x = ncread(full_dir, ['/', this_flowline, '/x']);
            y = ncread(full_dir, ['/', this_flowline, '/y']);
            nc_data.(this_glacier).(this_flowline) = [x,y];
        catch
            % if no data, then skip to the next
            continue
        end
    end
        
end
end

