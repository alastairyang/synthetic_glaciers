% Read in NetCDF files and output timeseries
clc
clear

%% Get all interested .nc files
% get the path to the uppermost directory for later scanning
original_dir = pwd;
objective_dir = 'ISMIP6 Data and processing';
% trim the path till right before objective_dir
position = strfind(pwd, objective_dir);
starting_path = original_dir(1:position+numel(objective_dir));
% back to where it was originally
cd(original_dir);

% we want both GIS (greenland ice sheet) and historical runs
key1 = "lithk"; % this should be the variable of interest
key2 = "GIS";
key3 = "historical";
key4 = ".nc";
key_remove = "MUN";
var = key1;
% retrieve an ambigous query keyword
keyword = combine_keys(key1, key2, key3, key4);

% query data from files
file_paths = search_files(starting_path, keyword);
file_new_paths = remove_paths(key_remove, file_paths);
ds_all = read_nc_var(file_new_paths, var);
cd(original_dir);

%% Get all interested .kml files
%%% import kml files which contain coordinates for all glaciers %%%
% The directory to the folder containing .kml data
kml_folder_path = '../../../Ice_velocity/kml_data';
% read all .kml files 
[kml_output, kml_names] = kml_data_out(kml_folder_path);
% convert to a table
glacier_info_table = struct2table(kml_output{1});

%% Get and print a list of experiments and indices
list_table = exp_list(ds_all) %#ok<NOPTS>

%% Read in and visualize flow line datasets
% The data coordinates are cartesian, hence no need for conversion
% for demonstration, let us first query NW
objective_dir = 'Data and SI';
% trim the path till right before objective_dir
position = strfind(pwd, objective_dir);
starting_path = original_dir(1:position+numel(objective_dir));
folder_name = 'flowlines_kristin/';
starting_path_folderadd = [starting_path, folder_name];

% search keywords
key1 = "WGreenland";
key2 = "ELA";
key3 = ".csv";
% unwanted keyword
key_remove = "Usurf";
keyword_csv = combine_keys(key1, key2, key3);

% read .csv flowline files
file_paths_csv = search_files(starting_path_folderadd, keyword_csv);
% remove entries with unwanted keyword
file_new_paths_csv = remove_paths(key_remove, file_paths_csv);

% get a cell array containing every flow line data in structures
% data sampling happens here
fl_out = flcsv2struct(file_new_paths_csv);

% for a given point (in cartesian), find the nearest five glaciers
% based off the data in fl_out
% data_interest is the point of interest
data_interest = [-155000, -2275000];
glacier_number = 5;
ng_out = nearest_glaciers(fl_out, data_interest, glacier_number);

%'VariableNames', {"X", "Y"}
% visualize data points
index = 10;
visualize_flowline(ds_all, ng_out, index);

cd(original_dir);
%% get time series
% convert lats and longs to polar stereographic
%glacier_psn = ll2psn(glacier_info_table.Y, glacier_info_table.X);
ts = point2timeseries(ng_out, ds_all, 1);

%% visualize
cellarray2plots(ts)

%%
%%%%% visualize greenland
pcolorpsn(ds_all{1,1}.lat, ds_all{1,1}.lon, reshape(ds_all{1,1}.var(:,:,1), size(ds_all{1,1}.lat)))
