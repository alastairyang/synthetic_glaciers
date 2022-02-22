function [kml_output, kml_names] = kml_data_out(folder_dir)
%   this file takes in a directory of a folder where all .kml files are
%   stored
%   Outputs are what are cells with data stored as struct 
    
    folder = folder_dir;
    if isa(folder, 'string')
        folder = convertStringToChars(folder);
    end
    
    % if the last char ends in "/"
    switch strcmp(folder(end),'/')
        case 1
            keyword = '*.kml';
        case 0
            keyword = '/*.kml';
    end
    
    path = [folder, keyword];
    file_paths = dir(path);
    % allocate space
    kml_output = cell(numel(file_paths), 1);
    kml_names = [];

    % read
    for i = 1:numel(file_paths)
        kml_output{i} = kml_shapefile([file_paths(i).folder,'/', file_paths(i).name]);
        sprintf('#%d .kml file imported!', i)
        
        if isa(file_paths(i).name, 'char')
            temp_name = convertCharsToStrings(file_paths(i).name);
        else
            temp_name = file_paths(i).name;
        end
        
        kml_names = [kml_names, temp_name];
    end
    
    
end

