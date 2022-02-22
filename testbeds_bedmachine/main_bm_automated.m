%% processing Felikson 2020 flowline data
clear 
clc

original_dir = pwd;
file_path = get_absolute_path('Research','clean_flowline_lengths.txt');
cd(original_dir)
ncfile_path = get_absolute_path('Research','flowline_Felikson_2020/netcdfs/');
cd(original_dir)

fl_idx_filedir = file_path{1};
fl_ncfolder_filedir = ncfile_path{1}(1:end-1);
fl_idx = flowline_fileindices(fl_idx_filedir);

% get all flow line coordinate data into nested structures
[glacier_idx_parsed, fl_data] = read_felikson_nc(fl_idx, fl_ncfolder_filedir);

%% Alternative, processing my own Thalweg data
clear
clc

original_dir = pwd;
file_path = get_absolute_path('Research','Thalweg_shapefile.shp');
cd(original_dir)

% read in the data
thalweg = read_thalweg_shp(file_path{1});


%% Get BedMachine data
original_dir = pwd;
objective_dir = 'Research';
% trim the path till right before objective_dir
position = strfind(pwd, objective_dir);
starting_path = original_dir(1:position+numel(objective_dir));
data_filename   = 'BedMachineGreenland-2021-04-20.nc';

% search for the bedmachine file
file_paths = search_files(starting_path, data_filename);
file_path = file_paths{1};

% now with the path, get bedmachine data
bm_data = import_bm(file_path);
X = bm_data.X;
Y = bm_data.Y;
x = bm_data.x;
y = bm_data.y;
b = bm_data.b;
h = bm_data.h;
s = bm_data.s;
mk = bm_data.mk;

cd(original_dir)

%% Extract geometry
method = 'straight';
map_data.X = X;
map_data.Y = Y;
map_data.b = b;
map_data.h = h;
map_data.s = s;
map_data.mk = mk;

% Visualize every flow line on a plot
%visualize_flowline_bm(fl_data)

% total 6 main(sub) flowlines per flowline
% 1: first, 3: the third, 6: the last
data_res = 'low';

if strcmp(data_res, 'low')
    data_all = flowline2geometry(thalweg, map_data, data_res, 1, method);
else % Dennis' data
    data_all = flowline2geometry(fl_data, map_data, data_res, 1, method);
end
%[data_crop, data_trans] = flowline2geometry(fl_data, map_data, 3, method);
%[data_crop, data_trans] = flowline2geometry(fl_data, map_data, 6, method);