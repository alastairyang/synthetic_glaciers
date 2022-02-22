function file_paths = search_files(search_dir, keyword)
%SEARCH evaluate if the current directory has the needed files
%   If there are, return the filenames and the relative path
%   If there are not, return the vector of all folder names
%   If there aren't even folders, print error
    

    % search dir should be an absolute path; it is the uppermost directory
    cd(search_dir)
    
    if isa(keyword, 'char')
        keyword = convertCharsToStrings(keyword);
    end
    % expected keyword form: **/keyword/ which can search in all subdirs
    full_keyword = "**/" + keyword;
    
    % keyword should be ambigous keyword recognisable by dir command
    all_files_struct = dir(full_keyword);
    
    N = numel(all_files_struct);
    all_namepath = cell(N, 1);
    for i = 1:numel(all_files_struct)
        % stitch together the folder dir and the file's name
        all_namepath{i} = fullfile(all_files_struct(i).folder, '/', all_files_struct(i).name);
    end
    
    file_paths = all_namepath;
end

