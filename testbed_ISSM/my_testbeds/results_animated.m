%% This script animates the glacier dynamic thickness change from multiple
% forcing experiments
%% Parameters
index = 0;
model_type  = 't';


%% Load the data
model_index = ['syn_', num2str(index)];
mr = load(['results/',model_index, '_', 'meltrates.mat']);
fr = load(['results/',model_index, '_', 'fric.mat']);
rB = load(['results/',model_index, '_', 'rheoB.mat']);
mr_fr = load(['results/', model_index, '_', 'meltrates_fric.mat']);
mr_rB = load(['results/', model_index, '_', 'meltrates_rheoB.mat']);
fr_rB = load(['results/',model_index, '_', 'rheoB_fric.mat']);
mr_fr_rB = load(['results/',model_index, '_', 'meltrates_rheoB_fric.mat']);

%% Main Script
% Both the surface elevation and elevation time series are animated
% synchronously
% we need to interp the meshgrid data back to our regular grid

gifname = 'meltrate_fric_rheoB.gif';
md = mr_fr_rB.md;

%% 
[geometry, ~] = query_data(index, model_type);
syn = testbed_data(geometry{1});
X = syn.X;
Y = syn.Y;
x = X(1,:);
y = Y(:,1);

nt = md.timestepping.final_time/md.timestepping.time_step;
t_selected = 1:floor(nt*0.05):nt;
if t_selected(end) ~= nt % make sure that the last time is always present
    t_selected = [t_selected, nt];
end
real_t_selected = 0.1*t_selected; % corresponding real year

% regrid and save the data
% Surface elevation, friction coefficient, rheology B, melt rates
elev = cell(size(t_selected));
melt = cell(size(t_selected));
rheo = cell(size(t_selected));
fric = cell(size(t_selected));

% check: if the size of the forcing data is the same as the spatial
% coordiate size ->  not time dependent
% hence we mannually a time slice in the middle of the sequence to the dataset
% Next step is to time-interpolate the spatial data (forcings)
if size(md.basalforcings.floatingice_melting_rate,1) == size(md.mesh.x,1)
    % constant forcing. Add t=0 and t=final sim time
    melt_ori = md.basalforcings.floatingice_melting_rate;
    melt_yr  = real_t_selected(floor(numel(real_t_selected)/2));
else
    melt_ori = md.basalforcings.floatingice_melting_rate(1:end-1,:);
    melt_yr  = md.basalforcings.floatingice_melting_rate(end,:);
end
if size(md.materials.rheology_B,1) == size(md.mesh.x,1)
    % add t=0 and t=final sim time
    rheo_ori = md.materials.rheology_B;
    rheo_yr  = real_t_selected(floor(numel(real_t_selected)/2));
else
    rheo_ori = md.materials.rheology_B(1:end-1,:);
    rheo_yr  = md.materials.rheology_B(end,:);
end
if size(md.friction.coefficient,1) == size(md.mesh.x,1)
    % add t=0 and t=final sim time
    fric_ori = md.friction.coefficient;
    fric_yr  = real_t_selected(floor(numel(real_t_selected)/2));
else
    fric_ori = md.friction.coefficient(1:end-1,:);
    fric_yr  = md.friction.coefficient(end,:);
end

% for easier interpolation (avoid using extrapolation), we check if the
% last year in the forcing data correspond to the last simulation year
% if no, we replicate the forcing from the last known year and add to a new
% data column (last sim year, hence constant extrapolation)
if melt_yr(end) ~= real_t_selected(end) ||...
   rheo_yr(end) ~= real_t_selected(end) ||...
   fric_yr(end) ~= real_t_selected(end)
    melt_ori = [melt_ori, melt_ori(:,end)];
    rheo_ori = [rheo_ori, rheo_ori(:,end)];
    fric_ori = [fric_ori, fric_ori(:,end)];
    % add the last year to year vectors
    melt_yr = [melt_yr, real_t_selected(end)];
    rheo_yr = [rheo_yr, real_t_selected(end)];
    fric_yr = [fric_yr, real_t_selected(end)];
end

% We do the same thing for the t(1)
if melt_yr(1) ~= real_t_selected(1) ||...
   rheo_yr(1) ~= real_t_selected(1) ||...
   fric_yr(1) ~= real_t_selected(1)
    melt_ori = [melt_ori(:,1), melt_ori];
    rheo_ori = [rheo_ori(:,1), rheo_ori];
    fric_ori = [fric_ori(:,1), fric_ori];
    % add the last year to year vectors
    melt_yr = [real_t_selected(1), melt_yr];
    rheo_yr = [real_t_selected(1), rheo_yr];
    fric_yr = [real_t_selected(1), fric_yr];
end

N_rows = size(melt_ori,1);
all_melt = zeros(N_rows, numel(t_selected));
all_rheo = zeros(N_rows, numel(t_selected));
all_fric = zeros(N_rows, numel(t_selected));
% Interpolate for each point
for i = 1:N_rows
    all_melt(i,:) = interp1(melt_yr, melt_ori(i,:), real_t_selected);
    all_rheo(i,:) = interp1(rheo_yr, rheo_ori(i,:), real_t_selected);
    all_fric(i,:) = interp1(fric_yr, fric_ori(i,:), real_t_selected);
end

disp('Spatial data time interpolation is complete')
% Now for each column (slice of time), we need to regrid data

% Sample points along the thalweg and derive time series
if rem(size(X,1), 2) == 0
    mid_i = size(X,1)/2;
else
    mid_i = (size(X,1)+1)/2;
end
sample_i = 1:5:100;
% find their coordinates
sample_x = X(1,sample_i);
sample_y = repmat(Y(mid_i,1), 1, numel(sample_i));

count_i = 0;
thalweg_sample_ht = [];
for i = t_selected
    count_i = count_i + 1;
    elev{count_i} = griddata(md.mesh.x, md.mesh.y, md.results.TransientSolution(i).Surface, X, Y);
    % we sequentially determine max and min of elev of all time, which we use for
    % colarbar limit later
    if count_i == 1
        elev_max = -1000; % small negative number
        elev_min = 1000; % big positive number
    else
        this_min = min(elev{count_i} - elev{count_i-1},[],'all');
        this_max = max(elev{count_i} - elev{count_i-1},[],'all');
    end
    if count_i ~= 1 && this_max > elev_max
        elev_max = this_max;
    end
    if count_i ~= 1 && this_min < elev_min
        elev_min = this_min;
    end
    thalweg_sample_h = elev{count_i}(mid_i, sample_i);
    thalweg_sample_ht = [thalweg_sample_ht; thalweg_sample_h];

    % Forcings
    % re-grid each forcing data
    melt{count_i} = griddata(md.mesh.x, md.mesh.y, all_melt(:,count_i), X, Y);
    fric{count_i} = griddata(md.mesh.x, md.mesh.y, all_fric(:,count_i), X, Y);
    rheo{count_i} = griddata(md.mesh.x, md.mesh.y, all_rheo(:,count_i), X, Y);
end

thalweg_sample_ht = thalweg_sample_ht';
thalweg_sample_ht = thalweg_sample_ht./thalweg_sample_ht(:,1);

disp('Regridding of forcing and elevation data is complete')


%%
pausetime = 0.2;
figure('Position',[100, 100, 1000,800]);

% Create color gradient
color_length = numel(sample_i);
red = [255, 51, 153]/255;
sth = [153, 153, 255]/255;
colors_p = [linspace(red(1),sth(1),color_length)',...
            linspace(red(2),sth(2),color_length)',...
            linspace(red(3),sth(3),color_length)'];

% start half way, skip the first 50 years
N_t_selected = numel(t_selected);
start = floor(N_t_selected/2);

% Plot and make it a .gif file
gif(gifname,'DelayTime',0.3)
count_i = 0;
for i = start:N_t_selected
    count_i = count_i + 1;
    % update title 
    titlestr = ['Year = ', num2str(real_t_selected(i))];
    sgtitle(titlestr)

    subplot(3,2,1)
    imagesc(x, y, fric{i})
    colorbar
    caxis([min(all_fric,[],'all'), max(all_fric,[],'all')])
    title('Frictional Coefficient')

    subplot(3,2,2)
    imagesc(x, y, rheo{i})
    colorbar
    caxis([min(all_rheo,[],'all'), max(all_rheo,[],'all')])
    title('Rheology B')

    subplot(3,2,3)
    imagesc(x, y, melt{i})
    colorbar
    caxis([min(all_melt,[],'all'), max(all_melt,[],'all')])
    title('Floating ice basal melting rates')

    % In this subplot, we plot the 2D elevation change map
    subplot(3,2,4)
    if count_i == 1
        imagesc(x, y, zeros(size(elev{i})))
    else
        imagesc(x, y, elev{i} - elev{i-1})
    end
    colorbar
    caxis([elev_min, elev_max])
    title('Elevation change dh/dt')
    hold on
    scatter(sample_x, sample_y, 20, colors_p,'filled')
    
    % Surface elevation time series
    subplot(3,2,5)
    colororder(colors_p)
    for j = 1:numel(sample_i)
        plot(real_t_selected(1:i), thalweg_sample_ht(j,1:i),...
             '-*'); hold on;
    end
    xlim([real_t_selected(start), real_t_selected(end)])
    ylim([0, 1.5])
    title('Normalized elevation change at sampled locations')

    pause(pausetime)
    gif
end

%% Check if forcings are independent
% figure; 
% subplot(1,2,1); colororder(colors_p); plot(real_t_selected, ht.mr_fr_rB,'-*'); ylim([0,1]); xlim([50,100]); 
% title('Simulated with all forcings on')
% xlabel('Years')
% ylabel('Normalized elevation change')
% subplot(1,2,2); colororder(colors_p); plot(real_t_selected, fr_rBnmr,'-*'); ylim([0,1]);xlim([50,100]);
% title('Individual forcing added together')
% xlabel('Years')
% ylabel('Normalized elevation change')