% debug file
% this file is contingent on running netCDF_to_Timeseries
%% verify the transformed polar stereographic is consistent with given X,Y
% file of interest, index in the file_paths cell array 
n = 16;
filepath = file_paths{n};

lon = ncread(filepath, 'lon');
lat = ncread(filepath, 'lat');
[X, Y] = ll2psn(lat, lon);

x_c = double(transpose(ncread(filepath, 'x')));
y_c = double(ncread(filepath, 'y'));
[X_c, Y_c] = meshgrid(x_c, y_c);
X_c = X_c';
Y_c = Y_c';

error_X = X_c - X;
error_Y = Y_c - Y;

%% Visualize greenland and added points with given coordinates
% use the GrIS coordinates from a given nc file

n = 16;
ds = ds_all{n};
% get transformed polar stereographic coordinates for the given data
glacier_psn = ll2psn(glacier_info_table.Y, glacier_info_table.X);
glacier_latlon = [glacier_info_table.Y, glacier_info_table.X];

% figure 1, with transformed polar stereographic using my own func
figure1 = figure;
% reshape value matrix to a vector for coloring
% the variable value is just one slice in time
var_plot = reshape(squeeze(ds.var(:,:,1)), size(ds.X));
s = pcolor(ds.X, ds.Y, var_plot);
s.FaceColor = 'interp';
s.EdgeColor = 'none';
hold on
scatter(glacier_psn(:,1), glacier_psn(:,2),'filled','r');
hold off

% % subplot2, with lon and lat.
% % REPEAT: reshape value matrix to a vector for coloring
% subplot(1,2,2)
% var_plot = reshape(ds.var(:,:,1), size(ds.X));
% pcolorpsn(ds.lat, ds.lon, var_plot)
% hold on
% scatterpsn(glacier_latlon(:,1), glacier_latlon(:,2),'filled','g')
% hold off
% axis equal

%%  Visualize the lon_bnds and lat_bnds in a nc file, if there is any
%LATLON_BOUNDS
n = 1;
filepath = file_paths{n};

    lat_bnds = ncread(filepath, 'lat_bnds');
    lon_bnds = ncread(filepath, 'lon_bnds');
    
    lat_bnds1 = squeeze(lat_bnds(1,:,:));
    lon_bnds1 = squeeze(lon_bnds(1,:,:));
    [X1, Y1] = ll2psn(lat_bnds1, lon_bnds1);
    lat_bnds2 = squeeze(lat_bnds(2,:,:));
    lon_bnds2 = squeeze(lon_bnds(2,:,:));
    [X2, Y2] = ll2psn(lat_bnds2, lon_bnds2);
    lat_bnds3 = squeeze(lat_bnds(3,:,:));
    lon_bnds3 = squeeze(lon_bnds(3,:,:));
    [X3, Y3] = ll2psn(lat_bnds3, lon_bnds3);
    lat_bnds4 = squeeze(lat_bnds(4,:,:));
    lon_bnds4 = squeeze(lon_bnds(4,:,:));
    [X4, Y4] = ll2psn(lat_bnds4, lon_bnds4);
    
    figure
    mesh(X1, Y1, ones(size(X1)))
    hold on
    mesh(X2, Y2, 5*ones(size(X2)))
    hold on
    mesh(X3, Y3, 10*ones(size(X3)))
    hold on
    mesh(X4, Y4, 15*ones(size(X4)))
    hold off

