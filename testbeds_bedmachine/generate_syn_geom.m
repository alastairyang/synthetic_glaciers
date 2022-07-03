%% This file generates synthetic glacier geometry
% The numeric values of certain geometric variables (e.g., fjord width) are
% from the "md_var_combination" file.

% This file generates the glacier geometry file in .mat format, which also
% contains the forcing perturbation data (better integration in progress),
% and an attribute table detailing all the parameters used in its
% generation.

%% Parameter
sens_run = 0; % no sensitivity run

%% Read in model combination data and create labels
mdvar_combs = readtable('/Users/donglaiyang/Library/Mobile Documents/com~apple~CloudDocs/Documents/Documents_Alastair’s_MacBook_Pro/Buffalo/Research/git_research/testbeds_bedmachine/md_var_combinations.csv');
N_var = size(mdvar_combs,2);
N_mds = size(mdvar_combs,1);
varnames = mdvar_combs.Properties.VariableNames;
vartypes = repelem(["int32"], length(varnames)); % repeat element
% create a new table for labels
mdvar_combs_label = table('Size', size(mdvar_combs),...
                          'VariableTypes', vartypes,...
                          'VariableNames', varnames);

% create dictionary for the variables and their label initials
initial_dict = create_initial_dict();
% get the sequence of the label initials
varname_labels = cellfun(@(key) initial_dict(key), varnames);

% creating labels
% lower value is 0, high values: just keep adding 1 to the previous one
% e.g., if only two values, [0,1]; three values, [0,1,2]
for i = 1:length(varnames)
    col_data = mdvar_combs.(varnames{i});
    col_data_uniq = unique(col_data);
    M = containers.Map(transpose(sort(col_data_uniq)), 0:(length(col_data_uniq)-1));
    col_data_idx = arrayfun(@(key) M(key), col_data);
    mdvar_combs_label.(varnames{i}) = col_data_idx;    
end

% combine index and initials to form labels
% these few lines are so elegant. I won't be able to come up on my own
mdvar_combs_label = num2str(mdvar_combs_label.Variables,'%d');
varname_labels_matrix = repmat(varname_labels, [N_mds,1]);
varname_labels_matrix = varname_labels_matrix';
mdvar_combs_label = mdvar_combs_label';
row_interleave = reshape([varname_labels_matrix(:) mdvar_combs_label(:)]',2*size(varname_labels_matrix,1), []);
labels = row_interleave(:)';
% now labels is a long string, we need to use regular expression to split
labels = regexp(labels, sprintf('\\w{1,%d}', N_var*2), 'match');

%% CREATE THE GEOMETRY !!!!!!!!!!!!!!
for i = 1:N_mds
    label = labels(i);
    label = label{:};
    var_table = mdvar_combs(i,:);
    var_table.('label') = label;
    fx_gen_syn_geom(var_table, sens_run);
end


%% FUNCTION:
function dict = create_initial_dict()
    keys = {'fjord_width', 'groundingline_depth', 'basalfric_law', 'background_friccoef'};
    vals = {'f',            'g',                 'l',             'b'};
    dict = containers.Map(keys, vals);
end

function fx_gen_syn_geom(var_table, sens_run)
    %% Unpack the input table    
    fjord_wid_single = var_table.('fjord_width');
    basalfric_law = var_table.('basalfric_law');
    gl_depth = var_table.('groundingline_depth');
    background_friccoef = var_table.('background_friccoef');
    label = var_table.('label');
    
    %% Main script
    % define parameters
    % Geometry parameters
    max_len = 400.0; % cutoff distance of flowline, in km
    n = 150; % sine wave sampling points
    ap_coef = 0.05; % filler wave amplitude
    elev_bench = gl_depth; % grounding line depth in meters, also zero benchmark line for parabola
    fjord_width = fjord_wid_single*2; % fjord valley width in meter. must be even numbers and multiples of 2n
    slidinglaw = basalfric_law; % sliding law. 0: hard bed. 1: soft bed.
    tot_wid = fjord_width*2; %  total width of this synthetic fjord
    tot_length = 100000.0;  % glacier length in meter
    shear_length = 30000.0; % length of shear margins in meters
    % fjordwater_length = 3000.0; % where ice is absent (did i use this?)
    bed_shearstress = 100000.0; % kPa
    background_basalfric = background_friccoef; % background basal friction coefficient, unit ?
    shelf_length = 2900.0; % length of ice shelf in meters
    front_height = 10.0; % height of ice front in meters

    % forcing parameters
    g = 9.81; % gravitational constant
    shear_B = 6.8067e7; % default value for B; later modified for shear margin weakening
    s_margin_width = 8; % length for each shear margin side, 3*150m
    s_margin_weakcoef = 0.3; % weaken to _% of the default value
    fric_coef_ampfactor = 0.3; % weaken to _% of the default fric coef
    slippatch_upstream_dist = 10000.0; % 10000 meter upstream of the grounding line
    slippatch_wid = 0.5*(fjord_width/2); % sigma in the gaussian function; half of fjord half width
    smb_constant = 0; % 0 m we, constant surface mass balance rate (when not using gradient method)

    % profile option: 'linear' or 'plastic'
    profile = 'plastic';

    % about the steep upward jump shape
    bottom_top_diff = 700; % 700 meter diff from the top to bottom
    bottom_top_long = 7500; % the distance along thalweg

    if tot_wid < fjord_width
        disp('Total width should be no smaller than fjord width')
        return;
    end
%     if rem(fjordwater_length,n) ~= 0
%         disp(' Change the fjordwater length to be a mutiples of spatial interval')
%         return
%     end

    %% Construct synthetic geometry with parabolic functions
    openfig('Plots/Try/mesh_fig_bh/glacier050_thalweg_xy_mesh_bh.fig');

    a = get(gca, 'Children');
    xdata = get(a, 'XData');
    ydata = get(a, 'YData');
    zdata = get(a, 'ZData');

    data.b = zdata{3};
    data.s = zdata{2};
    data.X = xdata{3};
    data.Y = ydata{3};

    %% Synethetic thalweg
    % get average of the middle 5 flow lines if odd, middle 4 if even
    % Here are some tunable parameters
    %   1. wavelength
    %   2. amplitude
    %   3. intervals between the overdeepening
    %   4. amplitudes for the said intervals (e.g., 1/20 of the average
    %   amplitudes of the overdeepening)
% 
%     b = data.b;
%     Ny = size(b,1);
%     s = data.s;
% 
%     % Extract some data along the centerlines
%     if rem(Ny, 2) == 0 % even number
%         % sample the middle 4 lines
%         b_sampled = b((Ny/2-2):(Ny/2+1),:);
%         b_sampled_mean = mean(b_sampled, 1);
%         s_sampled = s((Ny/2-2):(Ny/2+1),:);
%         s_sampled_mean = mean(s_sampled, 1);
%     else % sample the middle 3 lines
%         b_sampled = b(((Ny+1)/2-1):((Ny+1)/2+1),:);
%         b_sampled_mean = mean(b_sampled, 1);
%         s_sampled = s(((Ny+1)/2-1):((Ny+1)/2+1),:);
%         s_sampled_mean = mean(s_sampled, 1);
%     end
% 
%     x = data.X(1,:);

    %% STARTING BUILDING

    N_x_tot  = numel(0:n:tot_length);
    x = 0:n:(N_x_tot-1)*n;
    % concatenate
    thalweg = zeros(1, N_x_tot);

    % Smoothen the thalweg line
    sample_interval = 5;
    x_sample = x(1:sample_interval:numel(x));
    thalweg_sample = thalweg(1:sample_interval:numel(x));
    if x_sample(end) ~= x(end) % if last x is not already sampled in
        x_sample = [x_sample, x(end)]; % add last point for interpolation
        thalweg_sample = [thalweg_sample, thalweg(end)];
    end
    thalweg_interp = interp1(x_sample, thalweg_sample, x, 'spline');
    thalweg = thalweg_interp;

    if sum(isnan(thalweg),'all') > 0
        disp('There is NaN in thalweg, try changing the interpolation')
        return
    end

    % figure;
    % subplot(1,2,1)
    % plot(x, thalweg)
    % ylim([-600 1400])
    % subplot(1,2,2)
    % plot(x, thalweg_base)

    %% Synthetic fjord width
    % each element along the column is a function of
    %   1. distance to thalweg, y_hat; it is parabolic
    %   2. thalweg elevation, z
    %   3. distance to the terminus, x
    %       a. this is relevant if we want a variation of fjord width, e.g.,
    %       constriction

    % elev_bench: for parabola
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

    width_x = fjord_width*ones(numel(thalweg),1);

    % populate each column
    thalweg_i = (N_tot_wid+1)/2; % thalweg index location
    for i = 1:numel(thalweg)
        this_halfwid = width_x(i)/2;
        this_crosssec = syn_b(:,i);
        this_y = Y(:,i);
        % let two x-axis intersections be x1, x2
        % assuming x-axis centered parabolic curve
        % hence b = 0 (in ax^2+bx+c = 0), we find a = -(1/x1^2)*c
        % c = this_crosssec(thalweg_i)-elev_bench;
        c = this_crosssec(thalweg_i)-elev_bench;
        a = -(1/this_halfwid^2)*c;
        this_crosssec_formula = a.*this_y.^2 + c;
        syn_b(:,i) = this_crosssec_formula;
    end

    %% build fjord wall that levels off
    % make elements beyond the width zero
    wid_left_i  = thalweg_i - (fjord_width/n/2);
    wid_right_i = thalweg_i + (fjord_width/n/2);
    syn_b(1:wid_left_i-1,:) = 0;
    syn_b(wid_right_i+1:end,:)= 0;
    % check if it is an whole number
    if (floor(wid_left_i)~=wid_left_i) || (floor(wid_right_i)~=wid_right_i)
        disp('Index is not a whole number. Check if the width is an multiplier of 2n?')
        exit()
    end

    % find the numerical derivative at the edge of the overdeepening
    grad_od = zeros(numel(thalweg),1);
    for i = 1:numel(thalweg)
        grad_od(i) = (syn_b(wid_right_i,i)-syn_b(wid_right_i-2,i))/(2*n);
    end

    % modified error function: alpha*erf(beta*x)
    % its derivative at x=0 is 2*alpha*beta/sqrt(pi)
    % define alpha as the height of wall -> asymptotic value of the function
    alpha = 100+elev_bench; 
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

    %% Restore the elevation from elev_bench
    %syn_b = syn_b + elev_bench;

    %% Ice shelf base/ice thickness
    % need to data points to anchor a curve; bp stands for base point
    % density for flotation criterion
    rho_ice = 917; % kg/m^3
    rho_sw = 1023; % sea water, kg/m^3
    rho_hat = rho_ice/rho_sw;

    % provide a distance (x), then we find the base elevation there (on thalweg
    % line)
    bp_x = shelf_length;
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
    hp = [0, front_height]; % at the ice shelf front
    if hp(1) > ep(1)
        disp('the x of the head higher than the end now, choose a new x for head point')
        return
    end
    if hp(2) > ep(2)
        disp('head higher than end, choose a new end elevation value')
        return
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
    allnonzero_NaNbool = any(NaN_i, 2); % as long as this is element non-NaN
    gl_rows = find(allnonzero_NaNbool>0);
    % Search algorithm
    % start from the left end (ice front), find depth for each point until it
    % finds the first anchor point according to the flotation criterion
    gl_i = zeros(size(gl_rows));
    flot_pf = NaN(size(X));

    % This algorithm only works if the ice shelf TERMINATES AT THE END OF THE
    % DOMAIN (ice is present everywhere)
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
    % slope1 = 0.05;
    % slope2 = 0.025;
    % % the remaining x-length (after the ice shelf)
    % % divide in two parts
    % % the steeper part - 30% 
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
    % % spline interp the non-ice-shelf part
    % 
    % sample_interval = 20;
    % sample_i = 1:sample_interval:numel(x_syn_linear_noshelf); % n-point interval
    % x_syn_sample = x_syn_linear_noshelf(sample_i);
    % sl_syn_sample = sl_syn_linear(sample_i);
    % sl_syn_spline = spline(x_syn_sample, sl_syn_sample, x_syn_linear_noshelf);

    if strcmp(profile, 'plastic')
        % plastic rheology
        % H(x) = sqrt(2*ss*x/(rho*g))
        x_plastic = 0:n:(numel(x_syn_linear_noshelf)-1)*n;
        % displaced to the elevation at the end of the ice shelf
        sl_syn_spline = sqrt(2*bed_shearstress.*x_plastic/(rho_ice*g)) + s_end;

    elseif strcmp(profile, 'linear')
        % linear profile
        x_plastic = 0:n:(numel(x_syn_linear_noshelf)-1)*n;
        % we still use the end elevation from plastic profile to estimate
        % thickness near the ice divide
        s_plastic = sqrt(2*bed_shearstress.*x_plastic/(rho_ice*g)) + s_end;
        pf_slope = (s_plastic(end)-s_end)/(x_plastic(end)-x_plastic(1));
        sl_syn_spline = pf_slope*x_plastic + s_end;

    end

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

    %% Create ice shelf basal melting profile
    % % We approximate as a linear function of depth, where melt rate is higher
    % % around the grounding line.
    % % lower and upper bound of melt rate, meter/year;
    % % since it is Deep fjord Warm water type (DW), the melt rate at the bottom
    % % assumed to be higher than CR case, here set as 50 m/yr
    % range = [15, 50]; % melting rate range
    % flot_depth_range = [min(flot_pf,[],'all'), max(flot_pf,[],'all')];
    % unique_depth = unique(flot_pf);
    % unique_depth = sort(unique_depth(isfinite(unique_depth))); % increasing order
    % basal_melt = arrayfun(@(x) interp1(flot_depth_range, fliplr(range), x), flot_pf);
    % % let other NaN be 0, indicating no melt
    % basal_melt(isnan(basal_melt)) = 0;

    % Time dependent melting rate; time unit in years
    melt_years = {60, 80, 90}; % 0-1st, stationary; 1st-2nd, warming; 2nd-3rd, cooling
    melt_rates = [10, 90, 60];
    % range1 = range; % stationary period; specified already
    % range2 = [30, 100];
    % range3 = [20, 70];
    % 
    % basal_melt1 = basal_melt;
    % % 2nd
    % basal_melt2 = arrayfun(@(x) interp1(flot_depth_range, fliplr(range2), x), flot_pf);
    % basal_melt2(isnan(basal_melt2)) = 0;
    % % 3rd
    % basal_melt3 = arrayfun(@(x) interp1(flot_depth_range, fliplr(range3), x), flot_pf);
    % basal_melt3(isnan(basal_melt3)) = 0;

    basal_melt  = melt_rates(1)*ones(size(syn_b));
    basal_melt1 = basal_melt;
    basal_melt2 = melt_rates(2)*ones(size(syn_b));
    basal_melt3 = melt_rates(3)*ones(size(syn_b));

    % cell array, called transient_melt
    transient_melt = {basal_melt1, basal_melt2, basal_melt3};

    % Check consistency
    if numel(transient_melt) ~= numel(melt_years)
        disp('Melt year and rate are not consistent!')
        return;
    end

    % steady-state and spin-up melt rate
    constant_melt  = basal_melt;

    %% Friction coefficient map 1: constant and uniform map
    % % Here as a starter, we simply make it proportional to the topography,
    % % In the fjord trough, the friction is lower (also due to melt lubrication)
    % % than higher up in the wall
    % fric_coef_temp = syn_b;
    % range = [5, 200]; % lower and upper bound of coeff, as a reference to PIG example inverted range
    % normed_fric_coef = normalize(fric_coef_temp, 'range');
    % fric_coef = normed_fric_coef.*(range(2)-range(1));
    % fric_coef = fric_coef + range(1);
    fric_coef_cons = background_basalfric;
    fric_coef = fric_coef_cons*ones(size(syn_b));

    %% Friction coefficient map 2: time-dependent and slippery patch
    fric_coef_years = melt_years;
    fric_coef_ampfactors  = [1, fric_coef_ampfactor, 1]; % value between 0 and 1
    t_fric_coef   = cell(numel(fric_coef_years), 1);

    % specify the slippatch location
    % as the row of the thalweg, and the column as some distance upstream of the initial grounding line location.
    slippatch_yi = thalweg_i;
    slippatch_xi = find(isnan(flot_pf(thalweg_i,:)),1) + floor(slippatch_upstream_dist/150);
    slippatch_loc = [slippatch_xi, slippatch_yi];
    for i = 1:numel(fric_coef_years)
        % amplitude: if no perturbation, amp should be zero
        amp = fric_coef_ampfactors(i)*fric_coef_cons - fric_coef_cons;
        temp_fric_coef = transient_slippatch(fric_coef_cons, X, Y, slippatch_loc, slippatch_wid, amp);
        t_fric_coef{i} = temp_fric_coef;
    end
    % aggregate to structure
    transient_fric_coef.data  = t_fric_coef;
    transient_fric_coef.years = fric_coef_years;

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

    SMB_grad.smbref = smbref;
    SMB_grad.href   = href;
    SMB_grad.b_pos  = b_pos;
    SMB_grad.b_neg  = b_neg;

    % To avoid elevation-SMB feedback, we can simply use constant SMB rate
    SMB_cons = smb_constant*ones(size(syn_b));

    %% Shear margin 1: constant rheology B
    % first, we identify the location of the shear margin by looking at the
    % gradient of bed topography across y-axis
    syn_b_grad = transpose(gradient(syn_b'));
    % locating the shear margin: the row where max and min gradient values are
    % found
    [~, max_i] = max(syn_b_grad, [], 'all', 'linear');
    [~, min_i] = min(syn_b_grad, [], 'all', 'linear');
    [row_min, col_min] = ind2sub(size(syn_b_grad),min_i);
    [row_max, col_max] = ind2sub(size(syn_b_grad),max_i);

    % now get the index for both margins
    s_margin_i_left  = row_min:row_min + s_margin_width;
    s_margin_i_right = row_max:-1:row_max - s_margin_width;
    syn_rheoB = shear_B*ones(size(syn_b));
    syn_rheoB_unif = syn_rheoB;
    % substitute the weakening
    syn_rheoB_weak = syn_rheoB;
    rheoB_weak = shear_B*s_margin_weakcoef;
    syn_rheoB_weak(s_margin_i_left,1:floor(shear_length/150))  = rheoB_weak;
    syn_rheoB_weak(s_margin_i_right,1:floor(shear_length/150)) = rheoB_weak;
    % aggregate into a structure
    rheoB.rheoB_weak = syn_rheoB_weak;
    rheoB.rheoB_unif = syn_rheoB_unif;

    %% Shear margin 2: time dependent rheology B
    s_margin_years = melt_years;
    s_margin_amps  = [1, s_margin_weakcoef, 1]; % value between 0 and 1
    t_rheoB_weak   = cell(numel(melt_years), 1);
    for i = 1:numel(s_margin_years)
        temp_rheoB_weak = transient_shearmargin(shear_B, X, Y, s_margin_i_left, s_margin_i_right, s_margin_amps(i), shear_length);
        t_rheoB_weak{i} = temp_rheoB_weak;
    end
    % aggregate to structure
    transient_rheoB.data  = t_rheoB_weak;
    transient_rheoB.years = s_margin_years;

    % visualize the time-dependent rheoB map
    % figure; 
    % for i = 1:3
    %     imagesc(x,y,t_rheoB_weak{i}); 
    %     caxis([s_margin_amps(1)*shear_B, shear_B]);pause(2);
    % end
    %% Create multiple perturbed forcing datasets for sensitivity analysis
    if sens_run == 1
        N_sens = 3;
        % shelf basal melt rate
        max_melt_rates = linspace(max(melt_rates)-20, max(melt_rates)+20, N_sens);
        if any(max_melt_rates <= 0)
            disp(' Some melt rates have negative values')
            return
        end
        for i = 1:N_sens
            melt_rates = [melt_rates(1), max_melt_rates(i), melt_rates(3)];
            basal_melt1 = melt_rates(1)*ones(size(syn_b));
            basal_melt2 = melt_rates(2)*ones(size(syn_b));
            basal_melt3 = melt_rates(3)*ones(size(syn_b));

            % cell array, called transient_melt
            transient_melt = {basal_melt1, basal_melt2, basal_melt3};
            % add to the structure
            fieldname = ['test_',num2str(i)];
            sens_shelf_melt.(fieldname).transient_melt = transient_melt;
        end

        % frictional coefficient
        fric_coef_conss = linspace(fric_coef_cons-40, fric_coef_cons+40, N_sens);
        if any(fric_coef_conss <= 0)
            disp(' Some fric coefs have negative values')
            return
        end
        for i = 1:N_sens
            fieldname = ['test_', num2str(i)];
            sens_fric_coef.(fieldname) = fric_coef_conss(i)*ones(size(syn_b));
        end

        % rheology B, shear margin weakening
        s_margin_weakcoefs = linspace(s_margin_weakcoef, s_margin_weakcoef+0.4, N_sens);
        if any(s_margin_weakcoefs <= 0)
            disp(' Some rheoB have negative values')
            return
        end
        for i = 1:N_sens
            fieldname = ['test_', num2str(i)];
            sens_rheoB.(fieldname) = s_margin_weakcoefs(i)*shear_B*ones(size(syn_b));
        end
    end
%     %% visualize
%     figure('Position',[100, 100, 1500, 600])
%     subplot(2,3,1)
%     mesh(X,Y,syn_b)
%     title('Bed')
%     colorbar
% 
%     subplot(2,3,2)
%     mesh(X,Y,syn_s); hold on; mesh(X,Y,syn_base); hold off
%     title('Surface and Base')
%     colorbar
% 
%     subplot(2,3,3)
%     % plot lateral profile
%     plot(x, syn_base(thalweg_i, :),'r')
%     hold on
%     plot(x, syn_s(thalweg_i, :), 'b')
%     hold off
%     legend('Base','Surface')
% 
%     subplot(2,3,4)
%     imagesc(x,y,ocean_mask)
%     title('Floating/grounded Ice')
%     axis equal
%     xlim([min(x), max(x)])
%     ylim([min(y), max(y)])
%     colorbar
% 
%     subplot(2,3,5)
%     imagesc(x,y,transient_fric_coef.data{2});
%     title('Frictional Coef')
%     axis equal
%     xlim([min(x), max(x)])
%     ylim([min(y), max(y)])
%     colorbar
% 
%     subplot(2,3,6)
%     imagesc(x,y,transient_rheoB.data{2})
%     title('Rheology B Value')
%     axis equal
%     xlim([min(x), max(x)])
%     ylim([min(y), max(y)])
%     colorbar
% 
%     %% Only the geometry
%     figure;
%     mesh(X,Y,syn_s)
%     colormap summer; hold on
%     mesh(X,Y,syn_base); hold off
%     shading interp
%     grid off
%     axis off

    %% create an attribute table
    % building on top of the model variable table (the input of this
    % function)
    attrs_table = var_table;
    attrs_table.('fjord_width') = fjord_width; % update to the full width
    attrs_table.('ice_front_height') = front_height;
    attrs_table.('iceshelf_length') = shelf_length;
    attrs_table.('shearmargin_length') = shear_length;
    attrs_table.('shearmargin_width') = s_margin_width*150; 
    attrs_table.('shearmargin_weak_to_percent') = s_margin_weakcoef*100; % convert to percent
    attrs_table.('fric_coef_weak_to_percent') = fric_coef_ampfactor*100; % convert to percent
    attrs_table.('slippatch_upstream_dist') = slippatch_upstream_dist;
    attrs_table.('slippatch_width') = slippatch_wid;
    attrs_table.('parabolicprofile_bed_strength') = bed_shearstress;
    attrs_table.('total_length') = tot_length;
    attrs_table.('default_rheoB') = shear_B;
    attrs_table.('author') = 'Donglai Yang';
    attrs_table.('creation_date') = string(datetime('today'));
    %% Additional notes
    notes = 'sliding law represents the p number in Paterson sliding law. p=1 is linear; p=3 is nonlinear. Label 0 represents p=1; label 1 represents p=3';
    
    %% Save the geometry as a structure
    syn.X = X;
    syn.Y = Y;
    syn.bed = syn_b;
    syn.base = syn_base;
    syn.s = syn_s;
    syn.h = syn_h;
    syn.ocean_mask = ocean_mask;
    syn.ice_mask   = ice_mask;
    syn.end_year   = melt_years{end};
    syn.fric_coef  = fric_coef;
    syn.transient_fric_coef = transient_fric_coef;
    syn.shelf_melt.transient_melt = transient_melt;
    syn.shelf_melt.constant_melt  = constant_melt;
    syn.shelf_melt.melt_years     = melt_years;
    syn.SMB_grad = SMB_grad;
    syn.SMB_cons = SMB_cons;
    syn.rheoB = rheoB;
    syn.transient_rheoB = transient_rheoB;
    syn.slidingP = slidinglaw;
    syn.label = label;
    syn.attrs_table = attrs_table;
    syn.notes = notes;
    
    savepath = ['Synthetic glaciers/newlabels/geom_', label,'.mat'];
    save(savepath,'syn')
    disp(['Geometry ', label,' has been sucessfully created!! :) '])

    %% Save the perturbed forcings for sensitivity runs
    if sens_run == 1
        %sens_3.shelf_melt = sens_shelf_melt;
        sens.fric_coef  = sens_fric_coef;
        sens.rheoB      = sens_rheoB;
        % count the exisiting number of fields (must be before N_test)
        N_vars = numel(fieldnames(sens));
        sens.N_vars     = N_vars;

        sens.N_test     = N_sens;
        % save('Sensitivity/sens_0','sens_0')
    end
end