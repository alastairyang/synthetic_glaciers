%% This script visualize the timeseries (normalized h change)
% it plots one model results (e.g., all forcing active) from multiple
% synthetic glacier simulations and makes a 1xn figure
%% Parameters
figure('Position',[100,100,1120,200])
subplot(1,3,1)
index = '00';
model_type  = 't';
model_index = ['syn_', index];
all_00 = load(['results/',model_index, '/', 'meltrates_rheoB_fric.mat']);
md = all_00.md;
plot_ht(md, index, model_type)
xlabel('Time (year)','FontSize',13, 'FontName','Times')
ylabel('Elevation change wrt t=0','FontSize',13, 'FontName','Times')
title('Half width','FontSize',13, 'FontName','Times')

subplot(1,3,2)
index = '0';
model_type  = 't';
model_index = ['syn_', index];
all_0 = load(['results/',model_index, '/', 'meltrates_rheoB_fric.mat']);
md = all_0.md;
plot_ht(md, index, model_type)
xlabel('Time (year)','FontSize',13, 'FontName','Times')
title('Standard width','FontSize',13, 'FontName','Times')

subplot(1,3,3)
index = '01';
model_type  = 't';
model_index = ['syn_', index];
all_01 = load(['results/',model_index, '/', 'meltrates_rheoB_fric.mat']);
md = all_01.md;
plot_ht(md, index, model_type)
xlabel('Time (year)','FontSize',13, 'FontName','Times')
title('Double width','FontSize',13, 'FontName','Times')

print(gcf,'Graphs/timeseries.png','-dpng','-r300');  


%% APPENDIX: function
function [] = plot_ht(md, index, model_type)
    %% Main Script
    % Both the surface elevation and elevation time series are animated
    % synchronously
    % we need to interp the meshgrid data back to our regular grid
    % 
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
    base = cell(size(t_selected));
    Vel  = cell(size(t_selected));
    longi_dev_stress = cell(size(t_selected));
    later_dev_stress = cell(size(t_selected));

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
    sample_i = 1:10:150;
    % find their coordinates
    sample_x = X(1,sample_i);
    sample_y = repmat(Y(mid_i,1), 1, numel(sample_i));

    count_i = 0;
    thalweg_sample_ht = [];
    lateral_h = zeros(numel(t_selected), size(X,2));
    lateral_b = zeros(size(lateral_h));
    melt_t = zeros(size(t_selected));
    for i = t_selected
        count_i = count_i + 1;
        % surface elevation and base
        elev{count_i} = griddata(md.mesh.x, md.mesh.y, md.results.TransientSolution(i).Surface, X, Y);
        base{count_i} = griddata(md.mesh.x, md.mesh.y, md.results.TransientSolution(i).Base, X, Y);

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

        % Get the lateral profiles
        lateral_h(count_i,:) = elev{count_i}(mid_i, :);
        lateral_b(count_i,:) = base{count_i}(mid_i, :);

        % get the single melt rate value
        melt_t(count_i) = mean(melt{count_i},'all');

        % Get the velocity field Vel
        Vel{count_i} = griddata(md.mesh.x, md.mesh.y, md.results.TransientSolution(i).Vel, X, Y);
        % Sequentially determine velocity limits of all times
        this_vel_max = max(Vel{count_i},[],'all');
        this_vel_min = min(Vel{count_i},[],'all');
        if count_i == 1
            vel_max = this_vel_max;
            vel_min = this_vel_min;
        elseif this_vel_max > vel_max
            vel_max = this_vel_max;
        elseif this_vel_min < vel_min
            vel_min = this_vel_min;
        end

        % Get the longitudinal and lateral stress
        % first call 'mechanicalproperties' then grid data
        md = mechanicalproperties(md, md.results.TransientSolution(i).Vx,...
                                      md.results.TransientSolution(i).Vy);
        xx = md.results.deviatoricstress.xx;
        yy = md.results.deviatoricstress.yy;
        % Now xx and yy are defined at each element and covers the entire
        % element
        % We need to get the centroid locations and grid data
        [centroids_x, centroids_y] = triangle_centroid(md.mesh.elements, md.mesh.x, md.mesh.y);
        longi_dev_stress{count_i} = griddata(centroids_x, centroids_y, xx, X, Y);
        later_dev_stress{count_i} = griddata(centroids_x, centroids_y, yy, X, Y);   
        % we sequentially determine max and min of stress, which we use for
        % colarbar limit later
        this_longi_max = max(longi_dev_stress{count_i},[],'all');
        this_longi_min = min(longi_dev_stress{count_i},[],'all');
        this_later_max = max(later_dev_stress{count_i},[],'all');
        this_later_min = min(later_dev_stress{count_i},[],'all');
        if count_i == 1
            longi_max = this_longi_max;
            longi_min = this_longi_min;
            later_max = this_later_max;
            later_min = this_later_min;
            continue; % go to next iteration
        elseif this_longi_max > longi_max
            longi_max = this_longi_max;
        elseif this_longi_min < longi_min
            longi_min = this_longi_min;
        elseif this_later_max > later_max
            later_max = this_later_max;
        elseif this_later_min < later_min
            later_min = this_later_min;
        end

    end

    thalweg_sample_ht = thalweg_sample_ht';

    disp('Regridding of forcing and elevation data is complete')
    
    % Create color gradient
    color_length = numel(sample_i);
    red = [255, 51, 153]/255;
    sth = [153, 153, 255]/255;
    colors_p = [linspace(red(1),sth(1),color_length)',...
                linspace(red(2),sth(2),color_length)',...
                linspace(red(3),sth(3),color_length)'];

    % start half way, skip the first 50 years
    N_t_selected = numel(t_selected);
    start_i = floor(N_t_selected/2);

    % normalize surface elev time series wrt the starting point
    thalweg_sample_ht = thalweg_sample_ht./thalweg_sample_ht(:,start_i);
    
    % PLOT Surface elevation time series
    colororder(colors_p)
    count_i = 0;
    for i = start_i:N_t_selected
        count_i = count_i + 1;
        for j = 1:numel(sample_i)
            plot(real_t_selected(1:i), thalweg_sample_ht(j,1:i),...
                 '-*'); hold on;
        end
    end
    xlim([real_t_selected(start_i), real_t_selected(end)])
    ylim([0, 1.2])
end
