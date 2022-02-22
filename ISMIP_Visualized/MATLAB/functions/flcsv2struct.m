function fl_out = flcsv2struct(file_paths, n_sample)
%FLCSV2TENSORS Convert flow line csv files to a table with sampling due to
%a large number of data points
%
%   Example: fl_out = flcsv2struct(file_paths)
%
%   Input:
%       file_paths: cell arrays of file paths
%       n_sample: number of sampling points per flow line; if available
%                number of points is less, then incorporate all
%
%   Output:
%       fl_out: a table with X, Y, name
%       
    % total number of files to be processed
    N = numel(file_paths);
    % sampling target number
    N_target = n_sample;
    % Allocate memory
    coor = zeros(N*N_target, 2);
    name_v = strings(N*N_target, 1);
    for i = 1:N
        mx_temp = readmatrix(file_paths{i});
        
        % if there are more data points than N_target, we perform sampling
        if size(mx_temp, 1) > N_target
            rng(3);
            mx_temp = datasample(mx_temp, N_target, 'Replace', false);
        end
        
        % get the flow line name from the path/directiory
        % by locating the last left slash / and .csv
        full_dir = file_paths{i};
        sspositions = strfind(full_dir, "/");
        csvposition = strfind(full_dir, ".csv");
        % get from left slash onward
        name = full_dir(sspositions(end)+1 : csvposition-1);

        % number of points sampled: not always equal to N_target as some 
        % files have fewer points to begin with
        N_sampled = size(mx_temp, 1);
        % replicate the name for N_sampled times to make labels
        name_temp_v = repmat(convertCharsToStrings(name), N_sampled, 1);
        
        % append to the output matrix fl_out
        % find the starting and ending index for this segment of array
        idx_start = (i-1)*N_target + 1;
        idx_end   = idx_start + (N_sampled-1);
        
        coor(idx_start:idx_end, 1:2) = mx_temp;
        name_v(idx_start:idx_end, 1)   = name_temp_v;
    end
    % add name vector and make a table
    fl_out = table(coor(:,1), coor(:,2), name_v,...
                        'VariableNames', {'X','Y','Name'});
    
end

