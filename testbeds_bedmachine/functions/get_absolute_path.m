function file_paths = get_absolute_path(folder, file)
%GET_ABSOLUTE_PATH Get the absolute path to a file provided 
%   1. A folder name that exists in both pwd and its address
%   2. A file name
%
%   Example:
%       target address: /user/me/doc/file.txt
%       pwd:            /user/me/work/workfile.m
%       The target file is file.txt, and the folder name can be "me"

    original_dir = pwd;
    objective_dir = folder;
    % trim the path till right before objective_dir
    position = strfind(pwd, objective_dir);
    starting_path = original_dir(1:position+numel(objective_dir)-1);

    % keyword is the filename, normally
    keyword = file;
    
    % back to the original folder
    cd(original_dir)
    file_paths = search_files(starting_path, keyword);

end

