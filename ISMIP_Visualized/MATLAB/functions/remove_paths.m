function file_paths_new = remove_paths(key, file_paths)
%REMOVE_PATHS Summary remove entries from a file_paths cell array
%containing the given keywords
%   Input:
%       key: a string or char; keyword
%       file_paths: a cell array containing file paths for each cell
%   Output:
%       file_paths_new: new file paths cell array with entries removed

    if isa(key, 'char') % if it is char, it has only one keyword
        key = convertCharsToStrings(key);
    end
    if numel(key) > 1 % currently we only allow for one keyword
        error(' Currently we only allow for one keyword')
    end
    
    bool = ~contains(file_paths, key);
    file_paths_new = file_paths(bool);

end

