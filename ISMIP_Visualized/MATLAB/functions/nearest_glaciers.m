function ng_out = nearest_glaciers(g_data,interest_point, n)
%NEAREST_GLACIERS Find the coordinates of the nearest glaciers
%
%   Example: function ng_out = nearest_glaciers(g_data, interest_point, n)
%   
%   Input:
%       g_data: all glacier coordinates and labels in a table (X,Y,name)
%               expecting cartesian coordinates
%       interest_point: the coordinate of a point that is of interest to
%               you, also expecting cartesian coordinates
%       n: number of nearest glaciers you wish to sample
%   
%   Output:
%       ng_out: a table of selected glacier coordinates and labels
%       
    ng_points_out = [];
    ng_names_out  = [];
    i = 1;
    
    % first unpack table into a double array and string array
    coor  = table2array(g_data(:, 1:2));
    names = table2array(g_data(:, 3));
    
    while i <= n
        nearest_point = dsearchn(coor, interest_point);
        nearest_points_bool = ismember(names, names(nearest_point));
        nearest_points       = coor(nearest_points_bool, :);
        nearest_points_names = names(nearest_points_bool, :);
        % append
        ng_points_out = [ng_points_out; nearest_points]; 
        ng_names_out  = [ng_names_out; nearest_points_names];
        
        % now remove the ones that have been selected
        % meanwhile updating our original array
        coor  = coor(~nearest_points_bool, :);
        names = names(~nearest_points_bool, :);
        
        % report
        sprintf('Found No.%d glacier!', i)
        
        i = i+1;
        
    end
    
    % Combine two arrays into a table
    ng_out = table(ng_points_out(:,1), ng_points_out(:,2), ng_names_out, ...
                   'VariableNames', {'X', 'Y', 'Names'});
end

