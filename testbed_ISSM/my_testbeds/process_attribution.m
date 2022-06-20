%% Parameters
indices = ["00","0","01"]; % width gets bigger
N_idx = length(indices);
gcp_i = 1:10:150;
for i = 1:N_idx
    index = convertStringsToChars(indices(i));
    model_type  = 't';
    model_index = ['syn_', index];

    mr = load(['results/',model_index, '/', 'meltrates.mat']);
    md = mr.md;
    [hts.(['syn_',index]).mr, t] = get_ht(md, index, model_type, gcp_i);
    mr_diff = abs(1 - hts.(['syn_',index]).mr);

    fr = load(['results/',model_index, '/', 'fric.mat']);
    md = fr.md;
    [hts.(['syn_',index]).fr, ~] = get_ht(md, index, model_type, gcp_i);
    fr_diff = abs(1 - hts.(['syn_',index]).fr);

    rB = load(['results/',model_index, '/', 'rheoB.mat']);
    md = rB.md;
    [hts.(['syn_',index]).rB, ~] = get_ht(md, index, model_type, gcp_i);
    rB_diff = abs(1 - hts.(['syn_',index]).rB);

    mr_fr_rB = load(['results/',model_index, '/', 'meltrates_rheoB_fric.mat']);
    md = mr_fr_rB.md;
    [hts.(['syn_',index]).mr_fr_rB, ~] = get_ht(md, index, model_type, gcp_i);
    all_diff = abs(1 - hts.(['syn_',index]).mr_fr_rB);
    
    mr_pct = mr_diff./(mr_diff+fr_diff+rB_diff);
    mr_pct(:,t<60) = 1/3; % for the period w/o perturbation, set to 0.33
    fr_pct = fr_diff./(mr_diff+fr_diff+rB_diff);
    fr_pct(:,t<60) = 1/3;
    rB_pct = rB_diff./(mr_diff+fr_diff+rB_diff);
    rB_pct(:,t<60) = 1/3;
    
    % add to the field
    hts.(['syn_',index]).mr_pct = mr_pct;
    hts.(['syn_',index]).fr_pct = fr_pct;
    hts.(['syn_',index]).rB_pct = rB_pct;
    hts.(['syn_',index]).linear = linear;
    hts.(['syn_',index]).t      = t;
end

%% Make plots
samples_i = [1, 5, 10, 15];
N_samples = length(samples_i);

figure('Position',[100,100,1500,700]);
plot_count = 0;
plot_idx = [1:4;5:8;9:12;13:16];
geoplot_idx = plot_idx(:,1);
plot_idx = transpose(plot_idx(:,2:end));

for i = 1:N_samples
    sampled_i = samples_i(i);
    
    plot_count = plot_count+1;
    subplot(N_samples,4,plot_idx(plot_count))
    index = '00';
    plot(t, 100*hts.(['syn_',index]).mr_pct(sampled_i,:),'-*','LineWidth',3); hold on
    plot(t, 100*hts.(['syn_',index]).fr_pct(sampled_i,:),'-*','LineWidth',3); hold on
    plot(t, 100*hts.(['syn_',index]).rB_pct(sampled_i,:),'-*','LineWidth',3); hold off
    leg1 = legend('melt rate','fric. coef.','rheology B','Location','northwest');
    set(leg1,'Box','off')
    ylabel('Percentage','FontSize',13, 'FontName','Times')

    plot_count = plot_count+1;
    subplot(N_samples,4,plot_idx(plot_count))
    index = '0';
    plot(t, 100*hts.(['syn_',index]).mr_pct(sampled_i,:),'-*','LineWidth',3); hold on
    plot(t, 100*hts.(['syn_',index]).fr_pct(sampled_i,:),'-*','LineWidth',3); hold on
    plot(t, 100*hts.(['syn_',index]).rB_pct(sampled_i,:),'-*','LineWidth',3); hold off
    leg1 = legend('melt rate','fric. coef.','rheology B','Location','northwest');
    set(leg1,'Box','off')

    plot_count = plot_count+1;
    subplot(N_samples,4,plot_idx(plot_count))
    index = '01';
    plot(t, 100*hts.(['syn_',index]).mr_pct(sampled_i,:),'-*','LineWidth',3); hold on
    plot(t, 100*hts.(['syn_',index]).fr_pct(sampled_i,:),'-*','LineWidth',3); hold on
    plot(t, 100*hts.(['syn_',index]).rB_pct(sampled_i,:),'-*','LineWidth',3); hold off
    leg1 = legend('melt rate','fric. coef.','rheology B','Location','northwest');
    set(leg1,'Box','off')
end
% add title
subplot(N_samples,4,2)
title('Half width','FontSize',13, 'FontName','Times')
subplot(N_samples,4,3)
title('Standard width (7200 m)','FontSize',13, 'FontName','Times')
subplot(N_samples,4,4)
title('Double width','FontSize',13, 'FontName','Times')

% Add the plots showing the geometry and control points
% first query the geometry data
index = '0';
model_type = 't';
[geometry, ~] = query_data(index, model_type);
syn = testbed_data(geometry{1});
x = syn.X(1,:);
y = syn.Y(:,1);
h = syn.h;

% add to the three subplots
subplot(N_samples,4,geoplot_idx)
imagesc(y/1e3,x/1e3,h'); colormap(summer); hold on
xlabel('y (km)','FontSize',13, 'FontName','Times')
ylabel('x (km)','FontSize',13, 'FontName','Times')
a = colorbar;
a.Label.String = 'Ice Thickness (m)';
a.Label.FontSize = 13;
a.Label.FontName = 'Times';
% actually plotted control points
plotted_gcp_i = gcp_i(samples_i);
% color coding
color_length = length(plotted_gcp_i);
red = [255, 51, 153]/255;
sth = [153, 153, 255]/255;
colors_p = [linspace(red(1),sth(1),color_length)',...
            linspace(red(2),sth(2),color_length)',...
            linspace(red(3),sth(3),color_length)'];
        
if rem(size(syn.X,1), 2) == 0
    mid_i = size(syn.X,1)/2;
else
    mid_i = (size(syn.X,1)+1)/2;
end
sample_x = syn.X(1,plotted_gcp_i)/1e3;
sample_y = repmat(syn.Y(mid_i,1)/1e3, 1, numel(plotted_gcp_i));
scatter(sample_y, sample_x, 50, colors_p,'filled')

% Save to high resolution image
print(gcf,'Graphs/percentage.png','-dpng','-r300');  

% % add linearity plot
% subplot(N_samples+1,3,plot_count+1)
% index = '00';
% plot(t, hts.(['syn_',index]).linear(sampled_i,:),'-.','LineWidth',3); hold on
% 
% subplot(N_samples+1,3,plot_count+2)
% index = '0';
% plot(t, hts.(['syn_',index]).linear(sampled_i,:),'-.','LineWidth',3); hold on
% 
% subplot(N_samples+1,3,plot_count+3)
% index = '01';
% plot(t, hts.(['syn_',index]).linear(sampled_i,:),'-.','LineWidth',3); hold on

%% APPENDIX Functions
function [ht, t] = get_ht(md, index, model_type, sample_i)
    [geometry, ~] = query_data(index, model_type);
    syn = testbed_data(geometry{1});
    X = syn.X;
    Y = syn.Y;
    x = X(1,:);
    y = Y(:,1);

    nt = md.timestepping.final_time/md.timestepping.time_step;
    t_selected = 1:floor(nt*0.03):nt;
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

    % start half way, skip the first 50 years
    N_t_selected = numel(t_selected);
    start_i = floor(N_t_selected/2);

    % normalize surface elev time series wrt the starting point
    ht = thalweg_sample_ht./thalweg_sample_ht(:,start_i);
    
    ht = ht(:, start_i:end);
    t  = real_t_selected(start_i:end);

end
