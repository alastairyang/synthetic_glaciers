clc
clear

%% return to the 'Research' folder
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
mk = bm_data.mk;

cd(original_dir)

%% Experiment with the netCDF from Felikson et al 2020
foldername = 'flowline_Felikson_2020/netcdfs/';
testfile = 'glacier0001.nc';

gl1_fl3_d = ncread('glacier0001.nc', '/flowline03/d');
gl1_fl3_x = ncread('glacier0001.nc', '/flowline03/x');
gl1_fl3_y = ncread('glacier0001.nc', '/flowline03/y');

%% pick points
% finding grid limit

[xi, yi] = getline_zoom(h);
%[xi, yi] = getpts(fig1);

%% get mask
% re-assemble the map
map.X = X;
map.Y = Y;
map.b = b;
map.h = h;
vertices = [xi, yi];
% get the polygoned area (in a minimal bounding rectangle)
map_masked = polygon2mask(map, vertices);

% downsample for better visualization performance
[Xq, Yq, bq] = meshgrid_downsample(map_masked.X, map_masked.Y, map_masked.b);
[~,  ~,  hq] = meshgrid_downsample(map_masked.X, map_masked.Y, map_masked.h);
%% Rotation
theta = 40;
Rot = [cosd(theta) -sind(theta); sind(theta) cosd(theta)];
temp=[Xq(:),Yq(:)]*Rot' ;
      sz=size(Xq);
      
Xrot=reshape(temp(:,1),sz);
Yrot=reshape(temp(:,2),sz);

%% visualize
% plot GrIS ice thickness
plot_h()

%%
% only plot the masked region
figure;
mesh(Xrot, Yrot, bq)
view([30 60])
hold on
mesh(Xrot, Yrot, hq+bq)
hold off