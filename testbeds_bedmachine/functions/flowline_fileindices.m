function list = flowline_fileindices(filename)
%%This function extract the flow line indices for the corresponding NetCDF
%%files. Original data is from Felikson et al 2020.
%
%   Input: None
%
%   Output:
%       list: a cell array containing only the indices (no 'glacier'
%       prefix)
    
    T = readtable(filename);
    glacier_idx = T.(1);
    common_string = 'glacier';
    list = cellfun(@(x) x(strfind(x, common_string)+numel(common_string):end)...
                , glacier_idx, 'UniformOutput', false);
    
            
end