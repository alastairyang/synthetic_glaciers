%% define parameters
max_len = 30; % cutoff distance of flowline, in km
n = 50; % sine wave sampling points
ap_coef = 0.05; % filler wave amplitude
tot_wid = 4000; % 4 km total width of this synthetic fjord
elev_bench = 300; % zero benchmark line for parabola
width = 1200; % overdeepening width in meter. must be even numbers

%% Construct synthetic geometry with parabolic functions
openfig('Plots/mesh_fig_bh/glacier0003_flowline03_mesh_bh.fig');

a = get(gca, 'Children');
xdata = get(a, 'XData');
ydata = get(a, 'YData');
zdata = get(a, 'ZData');

% in zdata, {3} is base, {2} is 
data = re_resol_alongflow(zdata{3}, zdata{2}, max_len);

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
if rem(Ny, 2) == 0 % even
    bl1 = b(Ny/2-1,:);
    bl2 = b(Ny/2,  :);
    br1 = b(Ny/2+1,:);
    br2 = b(Ny/2+2,:);
    bav = (bl1+bl2+br1+br2)/4;
else
    bm  = b((Ny+1)/2,:);
    bl1 = b((Ny+1)/2-2,:);
    bl2 = b((Ny+1)/2-1,:);
    br1 = b((Ny+1)/2+1,:);
    br2 = b((Ny+1)/2+2,:);
    bav = (bm+bl1+bl2+br1+br2)/5;
end

% plot
%figure; mesh(data.X, data.Y, b);

% wavelength and amplitude from plot datapoints
wl0 = 4000;
wl1 = 17177.4-3145.16;
wl2 = 29758.1-24919.4;
ap0 = 150;
ap1 = -39.22-(-322.31);
ap2 = -57.82-(-355.02);
% od stands for overdeepening
% sample the sine wave; depend on wavelength
wave0_x = -wl0*0.75:50:0; % keep the latter 60% 
od0 = ap0*sin(pi/wl0.*wave0_x);
wave1_x = -wl1:50:0;
od1 = ap1*sin(pi/wl1.*wave1_x);
wave2_x = -wl2:50:0;
od2 = ap2*sin(pi/wl2.*wave2_x);

% low amplitude wave filling the intervals
wave_space1 = 2500;
ap_av = (ap0+ap1+ap2)/3;
% amplitude co-efficient
ap_sp = ap_coef*ap_av;
sp_wave_x = 0:50:wave_space1;
sp_wave   = ap_sp*sin(pi/wave_space1.*sp_wave_x);

% piecing together
N_x = numel(wave0_x) + numel(wave1_x) + numel(wave2_x) + 2*numel(sp_wave_x);
x = 0:50:(N_x-1)*n;
% concatenate
thalweg = [od0, sp_wave, od1, sp_wave, od2];

% Smoothen the thalweg line
x_sample = x(1:30:numel(x));
thalweg_sample = thalweg(1:30:numel(x));
thalweg_interp = interp1(x_sample, thalweg_sample, x, 'cubic');
thalweg = thalweg_interp;
% for od0, we don't want the curve down part; replace with straight line
[M,I] = min(thalweg(1:numel(od0)));
thalweg(1:I) = M;

if max(thalweg) > elev_bench % should adjust elev_bench
    disp('Should increase elev_bench value')
    exit
end

%% Synthetic fjord width
% each element along the column is a function of
%   1. distance to thalweg, y_hat; it is parabolic
%   2. thalweg elevation, z
%   3. distance to the terminus, x
%       a. this is relevant if we want a variation of fjord width, e.g.,
%       constriction

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
% create a plane of constant = this elevation bp_z
bp_plane = bp_z*ones(size(X));
bp_b_diff = bp_plane - syn_b;
bp_b_diff(bp_b_diff<0) = NaN;

%%%%% Specific to each glacier: find grounding line!! %%%%%
% This one, grounding line concave seaward

% first crop out the part beyond x = near_i 
% and find the rows for search
diff_crop = bp_b_diff(:,1:near_i);
NaN_bool = isfinite(diff_crop);
row_sum = sum(NaN_bool,2);
gl_rows = find(row_sum > 0);
% we now know the row numbers to start the search
gl_pos_i = zeros(size(gl_rows));

% search algorithm: starting at x = near_i
% find the first non-NaN
for j = 1:numel(gl_rows)
    row_num = gl_rows(j);
    this_row = diff_crop(row_num,:);
    success = 0;
    search_i = near_i;
    while success == 0
        if ~isnan(this_row(search_i))
            success = 1;
            gl_pos_i(j) = search_i;
        else
            search_i = search_i - 1;
        end
    end
end

% now with the positions obtained, get their exact bathymetric depth (due
% to meshgriding, they are not exactly the same)
surf = zeros(size(gl_pos_i));
depth = zeros(size(gl_pos_i));
thick = zeros(size(gl_pos_i));
for i = 1:numel(gl_pos_i)
    depth(i) = abs(syn_b(gl_rows(i), gl_pos_i(i)));
    thick(i) = depth(i)/rho_hat;
    surf(i) = thick(i) - depth(i);
end
% the shape of the ice shelf is that is linearly declines in surface elev.
min_gl_thick = min(surf,[],'all');
decline = (1-0.3)*min_gl_thick; % decline to 30 percent of surf elev
surf_elev_head = min_gl_thick - decline;

% again, iterate over all applicable rows to get the profiles
shelf_pf = NaN*ones(size(X));
for i = 1:numel(gl_pos_i)
    x1 = 0; 
    y1 = surf_elev_head;
    x2 = X(1,gl_pos_i(i)); 
    y2 = surf(i);
    xy = [x1,y1; x2,y2];
    x = X(1,1:gl_pos_i(i));
    temp = interp1(xy(:,1),xy(:,2),x);
    disp(numel(temp));
    shelf_pf(gl_rows(i),1:numel(temp)) = temp;
end

% Now we can get both ice shelf thickness profiles and base depth profiles
thick_all = shelf_pf.*(1/(1-rho_hat));
shelf_base = -(thick_all-shelf_pf);

% concatenate the surface elevation profile of non-marine ice
%% synthetic ice surface elevation; sl stnads for surface line, the profile line
s = data.s;
Ny = size(s,1);
slav = (s(70,:)+s(71,:)+s(72,:)+s(73,:))/4;

% we first iterate over all grouding line rows
syn_s = zeros(size(X));
for i = 1:numel(gl_pos_i)
    this_shelf = shelf_pf(gl_rows(i),:);
    p1 = [X(1,gl_pos_i(i)), this_shelf(gl_pos_i(i))];
    p2 = [5806.45, 450.00];
    p3 = [29758.1, 1038.82];
    p_syn = [p1; p2; p3];
    x_syn_linear_noshelf = X(1,gl_pos_i(i)+1:end);
    sl_syn_linear = interp1(p_syn(:,1),p_syn(:,2),x_syn_linear_noshelf);
    % spline interp the non-ice-shelf part
    
    sample_i = 1:50:numel(x_syn_linear_noshelf); % n-point interval
    x_syn_sample = x_syn_linear_noshelf(sample_i);
    sl_syn_sample = sl_syn_linear(sample_i);
    sl_syn_spline = spline(x_syn_sample, sl_syn_sample, x_syn_linear_noshelf);
    
    % concatenate together
    shelf_pf_finite = this_shelf(isfinite(this_shelf));
    syn_s(gl_rows(i), :) = [shelf_pf_finite, sl_syn_spline];
end

% Replicating the two marginal profiles to fill up the rest of the space
up_profile = syn_s(gl_rows(1),:);
bottom_profile = syn_s(gl_rows(end),:);
rep_up = repmat(up_profile, gl_rows(1)-1, 1);
rep_bottom = repmat(bottom_profile, size(X,1)-gl_rows(end), 1);
% do the same thing for shelf base
up_profile = shelf_base(gl_rows(1),:);
bottom_profile = shelf_base(gl_rows(end),:);
rep_up_b = repmat(up_profile, gl_rows(1)-1, 1);
rep_bottom_b = repmat(bottom_profile, size(X,1)-gl_rows(end), 1);
% substitute
syn_s(1:gl_rows(1)-1,:) = rep_up;
syn_s(gl_rows(end)+1:end,:) = rep_bottom;
% same for the shelf base
shelf_base(1:gl_rows(1)-1,:) = rep_up_b;
shelf_base(gl_rows(end)+1:end,:) = rep_bottom_b;
% NaN the extra part in shelf base
bs_b_diff = shelf_base - syn_b;
NaN_i = bs_b_diff<0;
shelf_base(NaN_i) = NaN;

% ice thickness = 0 -> 1 meter
s_b_diff = syn_s - syn_b;
syn_s(s_b_diff<0) = syn_b(s_b_diff<0) + 1;

figure; mesh(X,Y,syn_b);hold on;mesh(X,Y,syn_s)

%% Full synthetic ice thickness (floating tongue and anchored ice)
% replace the corresponding elements in base by base shelf
notNaN_i = isfinite(shelf_base);
b_syn_sub = syn_b;
b_syn_sub(notNaN_i) = shelf_base(notNaN_i);
syn_h = syn_s - b_syn_sub;

%% visualize
figure
mesh(X,Y,syn_b);hold on;mesh(X,Y,syn_s);hold on;mesh(X,Y,shelf_base)
%zlim([min(thalweg)-elev_bench, 0])

%% Save the geometry
syn_1.X = X;
syn_1.Y = Y;
syn_1.b = syn_b;
syn_1.s = syn_s;
syn_1.h = syn_h;
save('Synthetic glaciers/syn_1','syn_1')
