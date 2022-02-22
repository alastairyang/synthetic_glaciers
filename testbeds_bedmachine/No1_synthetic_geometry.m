%% define parameters
max_len = 400; % cutoff distance of flowline, in km
n = 150; % sine wave sampling points, same as BedMachine resolution
ap_coef = 0.05; % filler wave amplitude
tot_wid = 9000; % Arbitrary number; 4 km total width of this synthetic fjord
elev_bench = 300; % zero benchmark line for parabola
width = 4800; % overdeepening width in meter. must be even numbers
length = 45000; % total length of the fjord (along flow direction)
bed_shearstress = 100000; % 100kPa 
g = 9.81; % gravitational constant

if tot_wid < width
    disp('Total width should be no smaller than fjord width')
    return;
end

%% Construct synthetic geometry with parabolic functions
openfig('Plots/Try/mesh_fig_bh/glacier090_thalweg_xy_mesh_bh.fig');

a = get(gca, 'Children');
xdata = get(a, 'XData');
ydata = get(a, 'YData');
zdata = get(a, 'ZData');

data.b = zdata{3};
data.s = zdata{2};
data.X = xdata{3};
data.Y = ydata{3};

% % plot
% figure;
% mesh(data.X, data.Y, data.b); hold on; % bed
% mesh(data.X, data.Y, data.s); hold off

%% Synethetic thalweg
% get average of the middle 5 flow lines if odd, middle 4 if even
% Here are some tunable parameters
%   1. wavelength
%   2. amplitude
%   3. intervals between the overdeepening
%   4. amplitudes for the said intervals (e.g., 1/20 of the average
%   amplitudes of the overdeepening)

b = data.b;
Ny = size(b,1);
s = data.s;

% Extract some data along the centerlines
if rem(Ny, 2) == 0 % even number
    % sample the middle 4 lines
    b_sampled = b((Ny/2-2):(Ny/2+1),:);
    b_sampled_mean = mean(b_sampled, 1);
    s_sampled = s((Ny/2-2):(Ny/2+1),:);
    s_sampled_mean = mean(s_sampled, 1);
else % sample the middle 3 lines
    b_sampled = b(((Ny+1)/2-1):((Ny+1)/2+1),:);
    b_sampled_mean = mean(b_sampled, 1);
    s_sampled = s(((Ny+1)/2-1):((Ny+1)/2+1),:);
    s_sampled_mean = mean(s_sampled, 1);
end

x = data.X(1,:);
figure;
plot(x, b_sampled_mean,'r'); hold on;
plot(x, s_sampled_mean,'b'); hold off

%% START BUILDING THE GEOMETRY
% semi-wavelength and amplitude from plot datapoints
wl0 = 7000;
wl1 = 19244.2-6720.2;
wl2 = 36655.6-22298.8;
ap0 = 150;
ap1 = -39.22-(-322.31);
ap2 = -57.82-(-355.02);
% od stands for overdeepening
% sample the sine wave; depend on wavelength
wave0_x = -wl0:n:0; % keep the latter 60% 
od0 = ap0*sin(pi/wl0.*wave0_x);
wave1_x = -wl1:n:0;
od1 = ap1*sin(pi/wl1.*wave1_x);
wave2_x = -wl2:n:0;
od2 = ap2*sin(pi/wl2.*wave2_x);

% low amplitude wave filling the intervals
wave_space1 = 2500;
ap_av = (ap0+ap1+ap2)/3;
% amplitude co-efficient
ap_sp = ap_coef*ap_av;
sp_wave_x = 0:n:wave_space1;
sp_wave   = ap_sp*sin(pi/wave_space1.*sp_wave_x);

% piecing together
% concatenate
thalweg = [od0, sp_wave, od1, sp_wave, od2];

N_x_temp = numel(wave0_x) + numel(wave1_x) + numel(wave2_x) + 2*numel(sp_wave_x);
N_x_tot  = numel(0:n:length);
flat_filler = thalweg(end)*ones(1, N_x_tot - N_x_temp);
if N_x_tot < N_x_temp
    disp(' Overdeepening wavelength sum larger than total length!')
    return;
end
x = 0:n:(N_x_tot-1)*n;
% concatenate
thalweg = [od0, sp_wave, od1, sp_wave, od2, flat_filler];

% Smoothen the thalweg line
sample_interv = 30;
x_sample = x(1:sample_interv:numel(x));
thalweg_sample = thalweg(1:sample_interv:numel(x));
if x_sample(end) ~= x(end) % if last x is not already sampled in
    x_sample = [x_sample, x(end)]; % add last point for interpolation
    thalweg_sample = [thalweg_sample, thalweg(end)-1];
end
thalweg_interp = interp1(x_sample, thalweg_sample, x, 'cubic');
thalweg = thalweg_interp;
% for od0, we don't want the curve down part; replace with straight line
[M,I] = min(thalweg(1:numel(od0)));
thalweg(1:I) = M;

if sum(isnan(thalweg),'all') > 0
    disp('There is NaN in thalweg, try changing the interpolation')
    return
end

%% Synthetic fjord width
% each element along the column is a function of
%   1. distance to thalweg, y_hat; it is parabolic
%   2. thalweg elevation, z
%   3. distance to the terminus, x
%       a. this is relevant if we want a variation of fjord width, e.g.,
%       constriction


if max(thalweg) > elev_bench % should adjust elev_bench
    disp('Should increase elev_bench value')
    return;
end

% initialize a matrix where thalweg is the middle row
N_tot_wid = tot_wid/n+1; % odd number of rows
syn_b = zeros(N_tot_wid, numel(thalweg));
syn_b((N_tot_wid+1)/2,:) = thalweg;
% meshgrid coordinates for the synthetic
X = repmat(x, N_tot_wid, 1);
y = transpose(0:n:(N_tot_wid-1)*n) - (N_tot_wid-1)*n/2; % let y be symmetric wrt thalweg
Y = repmat(y, 1, numel(thalweg));

width_x = width*ones(numel(thalweg),1);

% populate each column
thalweg_i = (N_tot_wid+1)/2; % thalweg index location
for i = 1:numel(thalweg)
    this_halfwid = width_x(i)/2;
    this_crosssec = syn_b(:,i);
    this_y = Y(:,i);
    % let two x-axis intersections be x1, x2
    % assuming x-axis centered parabolic curve
    % hence b = 0, we find a = -(1/x1^2)*c
    c = this_crosssec(thalweg_i)-elev_bench;
    a = -(1/this_halfwid^2)*c;
    this_crosssec = a.*this_y.^2 + c;
    syn_b(:,i) = this_crosssec;
end

%% build fjord wall that levels off
% make elements beyond the width zero
wid_left_i  = thalweg_i - (width/n/2);
wid_right_i = thalweg_i + (width/n/2);
syn_b(1:wid_left_i-1,:) = 0;
syn_b(wid_right_i+1:end,:)= 0;

% find the numerical derivative at the edge of the overdeepening
grad_od = zeros(numel(thalweg),1);
for i = 1:numel(thalweg)
    grad_od(i) = (syn_b(wid_right_i,i)-syn_b(wid_right_i-2,i))/(2*n);
end

% modified error function: alpha*erf(beta*x)
% its derivative at x=0 is 2*alpha*beta/sqrt(pi)
% define alpha as the height of wall -> asymptotic value of the function
alpha = 100; % 100 m
beta = sqrt(pi)/(2*alpha).*grad_od;
N_side_wall = size(syn_b,1)-wid_right_i;
x_side_wall = transpose(1*n:n:N_side_wall*n);
syn_wall_right = zeros(numel(x_side_wall), numel(thalweg));
for i = 1:numel(thalweg)
    syn_wall_right(:,i) = alpha.*erf(beta(i).*x_side_wall);
end
syn_wall_left  = flipud(syn_wall_right);
% substitute the syn_b matrix
syn_b(1:wid_left_i-1,:) = syn_wall_left;
syn_b(wid_right_i+1:end,:) = syn_wall_right;

%% Ice shelf base/ice thickness
% need to data points to anchor a curve; bp stands for base point
% density for flotation criterion
rho_ice = 917; % kg/m^3
rho_sw = 1023; % sea water, kg/m^3
rho_hat = rho_ice/rho_sw;

% provide a distance (x), then we find the base elevation there (on thalweg
% line)
bp_x = 2800.0;
near_i = dsearchn(transpose(X(1,:)), bp_x);
bp_x = X(1,near_i);
bp_z = syn_b(thalweg_i, near_i);
bp = [bp_x,bp_z];
% find the required surface elevation there as a reference
s_bp = (-bp_z/rho_hat)-(-bp_z);
ep = [bp_x, s_bp]; % ep stands for end point, end point for shelf (potentially)
% create a synthetic ice shelf elevation that (temporarily) extends to the
% end of fjord length
% this profile is just a striaight line constrained by two points
hp = [0, 15]; % at the ice shelf front, the height is 15 m
if hp(2) > ep(2)
    disp('head higher than end, choose a new end elevation value')
    return;
end
% fjord end point, linear extrapolation
shelf_points = [hp;ep];
x_thalweg = X(1,:);
shelf_elev_line_temp = interp1(shelf_points(:,1), shelf_points(:,2), x_thalweg,...
                               'linear','extrap');
% make it a plane
shelf_elev_temp = repmat(shelf_elev_line_temp, size(X,1), 1);


%%%%% Specific to each glacier: find grounding line!! %%%%%
shelf_b_diff = shelf_elev_temp - syn_b;
shelf_elev_temp(shelf_b_diff<0) = NaN;
NaN_i = isfinite(shelf_elev_temp);
allnonzero_NaNbool = all(NaN_i, 2); % ice elevation must be always higher than base
gl_rows = find(allnonzero_NaNbool>0);
% Search algorithm
% start from the left end (ice front), find depth for each point until it
% finds the first anchor point according to the flotation criterion
gl_i = zeros(size(gl_rows));
flot_pf = NaN(size(X));

for j = 1:numel(gl_rows)
    row_num = gl_rows(j);
    this_shelfline = shelf_elev_temp(row_num,:);
    success = 0;
    search_i = 1;
    while success == 0
        % flotation criterion
        flot_depth = abs(this_shelfline(search_i)*(1/(1-rho_hat)-1));
        this_depth = abs(syn_b(row_num, search_i));
        if flot_depth < this_depth
            % record the depth for creating base profile later
            flot_pf(row_num, search_i) = -1*flot_depth;
            search_i = search_i + 1;
        else
            success = 1;
            % record the index and break the loop
            gl_i(j) = search_i;
            break;
        end
    end
end

%% Get the non-marine profile and replace the former temporay straight line pf.
% we start with the profile at the point of shelf profile most landward
[M,I] = max(gl_i); % most inland
x_end = X(gl_rows(I),M);
s_end = shelf_elev_temp(gl_rows(I),M);
this_shelfline = shelf_elev_temp(gl_rows(I),:);

% p1 = [x_end, s_end];
% slope1 = 0.034;
% slope2 = 0.028;
% % the remaining x-length (after the ice shelf)
% % divide in two parts
% x_rest = length - x_end;
% part1_portion = 0.3;
% x_part1 = x_rest*part1_portion;
% x_part2 = x_rest*(1-part1_portion);
% p2 = [x_end+x_part1, s_end+x_part1*slope1];
% p3 = [p2(1)+x_part2, p2(2)+x_part2*slope2];
% if p2(1)<x_end || p3(1)<p2(1) || p2(2)<s_end || p3(2)<p2(2)
%     disp(' Ice thickness profile is not ever-increasing.')
%     return
% end
% p_syn = [p1; p2; p3];
x_syn_linear_noshelf = X(1, M+1:end);
% sl_syn_linear = interp1(p_syn(:,1),p_syn(:,2),x_syn_linear_noshelf);
% spline interp the non-ice-shelf part
% 
% sample_i = 1:n:numel(x_syn_linear_noshelf); % n-point interval
% x_syn_sample = x_syn_linear_noshelf(sample_i);
% sl_syn_sample = sl_syn_linear(sample_i);
% sl_syn_spline = spline(x_syn_sample, sl_syn_sample, x_syn_linear_noshelf);

% plastic rheology
% H(x) = sqrt(2*ss*x/(rho*g))
x_plastic = 0:n:(numel(x_syn_linear_noshelf)-1)*n;
% displaced to the elevation at the end of the ice shelf
sl_syn_spline = sqrt(2*bed_shearstress.*x_plastic/(rho_ice*g)) + s_end;

% put together
this_shelfline(M+1:end) = sl_syn_spline;
syn_s = repmat(this_shelfline, size(X,1), 1);
s_b_diff = syn_s - syn_b;
syn_s(s_b_diff<0) = syn_b(s_b_diff<0) + 1;

%% synthetic ice surface elevation; sl stnads for surface line, the profile line
% ice thickness = 0 -> 1 meter
s_b_diff = syn_s - syn_b;
syn_s(s_b_diff<0) = syn_b(s_b_diff<0) + 1;

%% Ice thickness
shelf_base = flot_pf;
notNaN_i = isfinite(shelf_base);
b_syn_sub = syn_b;
b_syn_sub(notNaN_i) = shelf_base(notNaN_i);
syn_h = syn_s - b_syn_sub;

%% Glacier base
syn_base = flot_pf;
syn_base(isnan(flot_pf)) = syn_b(isnan(flot_pf));

%% Create mask
% Ocean_mask: Grounded vs Marine
ocean_mask = ones(size(shelf_base));
ocean_mask(notNaN_i) = 0; % negative: floating

% Ice_mask : ice presence or absence
% present: negative or -1
% absent:  positive or 1
ice_mask = -1.*ones(size(syn_h));
ice_mask(syn_h == 0) = 1;

%% Create friction coefficient map
% Here as a starter, we simply make it proportional to the topography,
% In the fjord trough, the friction is lower (also due to melt lubrication)
% than higher up in the wall
% fric_coef_temp = syn_b;
% range = [5, 200]; % lower and upper bound of coeff, as a reference to PIG example inverted range
% normed_fric_coef = normalize(fric_coef_temp, 'range');
% fric_coef = normed_fric_coef.*(range(2)-range(1));
% fric_coef = fric_coef + range(1);
fric_coef = 100*ones(size(syn_b));

%% Create ice shelf basal melting profile
% We approximate as a linear function of depth, where melt rate is higher
% around the grounding line.
% range = [15, 30]; % lower and upper bound of melt rate, meter/year
% flot_depth_range = [min(flot_pf,[],'all'), max(flot_pf,[],'all')];
% unique_depth = unique(flot_pf);
% unique_depth = sort(unique_depth(isfinite(unique_depth))); % increasing order
% basal_melt = arrayfun(@(x) interp1(flot_depth_range, fliplr(range), x), flot_pf);
% % let other NaN be 0, indicating no melt
% basal_melt(isnan(basal_melt)) = 0;
% 
% % Time dependent melting rate; time unit in years
 melt_years = {60, 80, 90}; % 0-1st, stationary; 1st-2nd, warming; 2nd-3rd, cooling
% range1 = range; % stationary period; specified already
% range2 = [30, 90];
% range3 = [20, 50];
% 
% basal_melt1 = basal_melt;
% % 2nd
% basal_melt2 = arrayfun(@(x) interp1(flot_depth_range, fliplr(range2), x), flot_pf);
% basal_melt2(isnan(basal_melt2)) = 0;
% % 3rd
% basal_melt3 = arrayfun(@(x) interp1(flot_depth_range, fliplr(range3), x), flot_pf);
% basal_melt3(isnan(basal_melt3)) = 0;

basal_melt  = 10*ones(size(syn_b));
basal_melt1 = basal_melt;
basal_melt2 = 60*ones(size(syn_b));
basal_melt3 = 45*ones(size(syn_b));

% cell array, called transient_melt
transient_melt = {basal_melt1, basal_melt2, basal_melt3};

% Check consistency
if numel(transient_melt) ~= numel(melt_years)
    disp('Melt year and rate are not consistent!')
    return;
end

% steady-state and spin-up melt rate
constant_melt  = basal_melt;


%% Create SMBgradient reference
% Elevation and SMB
% Data in reference to Helsen et al., 2012, 
% Coupling of climate models and ice sheet models by surface mass balance gradients: application to the Greenland Ice Sheet

% 1.href is just elevation or syn_s, unit is m
href = syn_s;
% 2.smbref is SMB calculated based on this referenece elevation field
% unit is kg m-2 yr-1, or equivalently, mm yr-1 water equivalent
ablation_p1 = [0, -2800]; 
ablation_p2 = [1400, 0];
ablation_elev = [ablation_p1(1), ablation_p2(1)];
ablation_SMB  = [ablation_p1(2), ablation_p2(2)];
accumulation_p1 = [1400, 0];
accumulation_p2 = [2000, 400];
accumulation_elev = [accumulation_p1(1), accumulation_p2(1)];
accumulation_SMB  = [accumulation_p1(2), accumulation_p2(2)];
ELA = accumulation_p1(1);
s_ablation = NaN*ones(size(syn_s));
index = find(syn_s < ELA);
s_ablation(index) = syn_s(index);
s_accumula = NaN*ones(size(syn_s));
index = find(syn_s >= ELA);
s_accumula(index) = syn_s(index);

% map elevation to their SMB
SMB_ablation = arrayfun(@(x) interp1(ablation_elev, ablation_SMB,x,...
                        'linear','extrap'), s_ablation);
SMB_accumula = arrayfun(@(x) interp1(accumulation_elev, accumulation_SMB, x,...
                        'linear','extrap'), s_accumula);
SMB_ablation(isnan(SMB_ablation)) = 0;
SMB_accumula(isnan(SMB_accumula)) = 0;
smbref = SMB_ablation + SMB_accumula;
                    
% obtain slopes for both regressions, to be used in ISSM
% b_pos is accumulation regime
temp = accumulation_p2 - accumulation_p1;
b_pos = temp(2)/temp(1);
% b_neg is ablation regime
temp = ablation_p2 - ablation_p1;
b_neg = temp(2)/temp(1);

SMB.smbref = smbref;
SMB.href   = href;
SMB.b_pos  = b_pos;
SMB.b_neg  = b_neg;

%% visualize
figure
subplot(2,1,1)
mesh(X,Y,syn_s); hold on; mesh(X,Y,syn_base);
hold on;
mesh(X,Y,syn_b); hold off
subplot(2,1,2)
mesh(data.X, data.Y, data.b); hold on; % bed
mesh(data.X, data.Y, data.s); hold off


%% Save the geometry
syn_1.X = X;
syn_1.Y = Y;
syn_1.b = syn_b;
syn_1.s = syn_s;
syn_1.h = syn_h;
syn_1.ocean_mask = ocean_mask;
syn_1.ice_mask   = ice_mask;
syn_1.end_year   = melt_years{end};
syn_1.fric_coef = fric_coef;
syn_1.shelf_melt.transient_melt = transient_melt;
syn_1.shelf_melt.constant_melt  = constant_melt;
syn_1.shelf_melt.melt_years     = melt_years;
syn_1.SMB = SMB;
syn_1.description = 'Calving on the Ridge type (CR) as described in Wood et al., 2021';
save('Synthetic glaciers/syn_1','syn_1')
