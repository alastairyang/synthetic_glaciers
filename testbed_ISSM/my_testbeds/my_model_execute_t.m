function output = my_model_execute_t(geometry_path, velocity_path, model_path, model_index, model_type)
%%  Model 
    tb = table();
    % the parameters() function is simply to store hyperparamters of
    % this model; they are not physical parameters
    params = parameters(tb);
    n_layer = params.n_layer;
    n_process = params.n_process;
    exponent = params.exponent;
    dt = params.dt;
    nt_t = params.nt_t;
    nt_spinup = params.nt_spinup;
    max_stress_grounded = params.max_stress_grounded;
    max_stress_floating = params.max_stress_floating;

%   If this is a spinup run, we set up model from scratch
    if strcmp(model_type, 'spinup')

        md = model;
        % parse out params
        % calling parameters(1) is to input something other than a table
        % so nothing is changed (changes should be made in the master file)
        % and hence this function is just to output

        % we use triangle method to create mesh
        % Thi is modeled after JI example
        % triangle(model, domain_file, average_element_size_meter)
        syn = testbed_data(geometry_path);
        [Xq,Yq,~] = meshgrid_downsample(syn.X, syn.Y, ones(size(syn.X))); % use ones as placeholder

        % generate Domain.exp file
        meshgrid2outline(Xq,Yq);
        % get a preliminary mesh which we will refine later
        md = bamg(md,'domain', 'Domain.exp', 'hmax', 3000);

        % Velocity field
        load(velocity_path); % loaded as 'V'
        vel_mesh = InterpFromGridToMesh(syn.x',syn.y,V.vel, md.mesh.x,md.mesh.y,0);
        % Speed
        md = bamg(md,'hmin',200,'hmax',2000,'field',vel_mesh,'err',5);
        %md=triangle(md, 'Domain.exp');
        plotmodel(md, 'data','mesh')

        % first, check what the vel looks like
        % plotmodel(md,'data',vel,'edgecolor','w')

    %Masks #2
        groundedice = double(InterpFromGridToMesh(syn.x',syn.y,syn.ocean_mask,md.mesh.x,md.mesh.y,0));
        groundedice(groundedice<=0) = -1; % when 0, make it negative, hence no ice
        icepresence=double(InterpFromGridToMesh(syn.x',syn.y,syn.ice_mask  ,md.mesh.x,md.mesh.y,0));

        %fill in the md.mask structure
        md.mask.ocean_levelset  = groundedice; %ice is grounded for mask equal one
        md.mask.ice_levelset    = icepresence;%ice is present when negatvie, and it is present everywhere
        % levelset method mandatory: levelset.spclevelset. Same convention
        % as ice_levelset


        %ploting
        plotmodel(md,'data','mesh', ...
                     'data',md.mask.ocean_levelset,'title','grounded/floating',...
                     'data',md.mask.ice_levelset,'title','ice/no-ice')

    %Parameterization #3
        ParamFile = ['parameters/syn_',num2str(model_index),'_',model_type, '.par'];
        md = parameterize(md, ParamFile);

        %Extrusion #?
        % only 5 layers exponent 1
        %md = extrude(md, n_layer, exponent);
        
        % rheology - B: data{1} should be uniform B, no shear margin
        % weakening. See it is using the first cell (no weakening)
        md.materials.rheology_B = double(InterpFromGridToMesh(syn.x',syn.y,syn.transient_rheoB.data{1},md.mesh.x,md.mesh.y,0));

        % friction coef
        % no perturbation (slippery patch)
        md.friction.coefficient = double(InterpFromGridToMesh(syn.x',syn.y,syn.transient_fric_coef.data{1},md.mesh.x,md.mesh.y,0));
        % no friction applied on floating ice; we did this before in .par
        % but a little redundancy does not hurt
        pos=find(md.mask.ocean_levelset<0);
        md.friction.coefficient(pos)=0;
        
        % forcing: needs to be forced after extrusion for model consistency
        % Ice shelf basal melt: constant melt rate
        % Constant melt rate forcing
        shelf_melt_mesh = InterpFromGridToMesh(syn.x', syn.y, syn.shelf_melt.transient_melt{1}, md.mesh.x, md.mesh.y, 0);
        md.basalforcings.floatingice_melting_rate = shelf_melt_mesh;

%         % adding time-dependent melt rate forcing
%         % THIS IS TEMPORARY
%         [x, y, X, Y, ~, ~, ~, ~, ~, ~, ~, ~, ~, shelf_melt, ~] = testbed_data(geometry_path);
%         N_years = numel(shelf_melt.melt_years);
%         for i = 1:N_years
%             shelf_melt_mesh = InterpFromGridToMesh(x', y, shelf_melt.transient_melt{i}, md.mesh.x, md.mesh.y, 0);
%             temp = [shelf_melt_mesh; shelf_melt.melt_years{i}]; % vertical concat year
%             if i == 1
%                 md.basalforcings.floatingice_melting_rate = temp;
%             else
%                 md.basalforcings.floatingice_melting_rate = [md.basalforcings.floatingice_melting_rate, temp];
%             end
%         end
        

%         % use PICO model for basal shelf melting rate
%         disp(' Start parameterizing shelf basal melting rate')
%         md.basalforcings = basalforcingspico();
%         md.basalforcings.basin_id = ones(md.mesh.numberofelements,1);
%         md.basalforcings.num_basins = 1;
%         md.basalforcings.maxboxcount = 5;
%         melt_year_axis = [1, 70, 90]; % year
%         md.basalforcings.farocean_temperature = [1.0 + 273.15, 10.0 + 273.15, 5.0 + 273.75; melt_year_axis];
%         md.basalforcings.farocean_salinity = [34.73, 34.73, 34.73; melt_year_axis]; % salinity; value taken from ISSM manual example as placeholder
%         md.basalforcings.isplume = 0;
%         md.basalforcings.geothermalflux = zeros(md.mesh.numberofvertices, 1);
%         md.basalforcings.groundedice_melting_rate = 0.0*ones(md.mesh.numberofvertices, 1);
        
        %set the surface mass balance
%         disp(' Start parameterizing SMB')
%         md.smb = SMBgradients();
%         href_mesh   = InterpFromGridToMesh(x', y, SMB_grad.href, md.mesh.x, md.mesh.y, 0);
%         smbref_mesh = InterpFromGridToMesh(x', y, SMB_grad.smbref, md.mesh.x, md.mesh.y, 0);
%         md.smb.href   = href_mesh;
%         md.smb.smbref = smbref_mesh;
%         % b_pos and b_neg also need to be fields
%         md.smb.b_pos = repmat(SMB_grad.b_pos, md.mesh.numberofvertices,1); 
%         md.smb.b_neg = repmat(SMB_grad.b_neg, md.mesh.numberofvertices,1);
        
        % if using constant smb
        %md.smb.mass_balance = InterpFromGridToMesh(x', y, SMB_cons, md.mesh.x, md.mesh.y, 0);
        
        % ADD HERE if other forcings are to be time-dependent
        % calving
        disp(' Start parameterizing calving')
        md.calving = calvingvonmises();
        md.calving.stress_threshold_groundedice = max_stress_grounded;
        md.calving.stress_threshold_floatingice = max_stress_floating;
        md.calving.min_thickness = 1;
        
        % frontal melting rate
        md.frontalforcings.meltingrate = 0.0*ones(md.mesh.numberofvertices, 1);
        
    %Set the flow computing method #5
        md = setflowequation(md, 'SSA', 'all');
        % before running the model, take a look at all boundary conditions
        plotmodel(md,'data','BC')

    %Solving #8
        md.cluster=generic('name',oshostname(),'np',n_process);
        %md.verbose=verbose('convergence',true);
        
%%%%%%%%%
    elseif strcmp(model_type, 't') % If this is a transient run with time-dependent climate forcing
        load(model_path)
        
        % This .par substitutes much of the final state of the spin-up run
        % to the initial conditions of this transient run
        ParamFile = ['parameters/syn_',num2str(model_index),'_',model_type, '.par'];
        md = parameterize(md, ParamFile);
        
        % adding time-dependent melt rate forcing
        syn = testbed_data(geometry_path);
        
        % time-dependent shear margin weakening
        N_years = numel(syn.transient_rheoB.years);
        for i = 1:N_years
            rheoB_weak_mesh = double(InterpFromGridToMesh(syn.x',syn.y,syn.transient_rheoB.data{i},md.mesh.x,md.mesh.y,0));
            temp = [rheoB_weak_mesh; syn.transient_rheoB.years{i}]; % vertical concat year
            if i == 1
                md.materials.rheology_B = temp;
            else
                md.materials.rheology_B = [md.materials.rheology_B, temp];
            end
        end

        % time-dependent frictional coefficient
        N_years = numel(syn.transient_fric_coef.years);
        for i = 1:N_years
            fric_coef_mesh = InterpFromGridToMesh(syn.x', syn.y, syn.transient_fric_coef.data{i}, md.mesh.x, md.mesh.y, 0);
            temp = [fric_coef_mesh; syn.transient_fric_coef.years{i}]; % vertical concat year
            if i == 1
                md.friction.coefficient = temp;
            else
                md.friction.coefficient = [md.friction.coefficient, temp];
            end
        end
        
        
        % time-dependent melt rate
        N_years = numel(syn.shelf_melt.melt_years);
        for i = 1:N_years
            shelf_melt_mesh = InterpFromGridToMesh(syn.x', syn.y, syn.shelf_melt.transient_melt{i}, md.mesh.x, md.mesh.y, 0);
            temp = [shelf_melt_mesh; syn.shelf_melt.melt_years{i}]; % vertical concat year
            if i == 1
                md.basalforcings.floatingice_melting_rate = temp;
            else
                md.basalforcings.floatingice_melting_rate = [md.basalforcings.floatingice_melting_rate, temp];
            end
        end      
        
        % Constant melt: just to make sure that this imported model stays
        % in steady-state
%         [x, y, X, Y, ~, ~, ~, ~, ~, ~, ocean_mask, ice_mask, end_year, ~, shelf_melt, SMB] = testbed_data(geometry_path);
%         shelf_melt_mesh = InterpFromGridToMesh(x', y, shelf_melt.transient_melt{1}, md.mesh.x, md.mesh.y, 0);
%         md.basalforcings.floatingice_melting_rate = shelf_melt_mesh;

%         md.basalforcings.groundedice_melting_rate = 2*ones(md.mesh.numberofvertices,1);

%         % PICO melt model
%         melt_year_axis = [10, 22, 28]; % year
%         md.basalforcings.farocean_temperature = [1.0 + 273.15, 6.0 + 273.15, 3 + 273.75; melt_year_axis];
%         md.basalforcings.farocean_salinity = [34.73, 34.73, 34.73; melt_year_axis]; % salinity; value taken from ISSM manual example as placeholder
        
        % Change mask numerical values (but not their classification)
        md.mask.ocean_levelset(md.mask.ocean_levelset<0) = -1;
        md.mask.ocean_levelset(md.mask.ocean_levelset>0) = 1;
        
    else % it is 't_sensitive' model type
        % For sensitivity analysis simulation, we are going to first
        % inherit the model from spin-up, then 
        load(model_path)
        
        % This .par substitutes much of the final state of the spin-up run
        % to the initial conditions of this transient run
        ParamFile = ['parameters/syn_',num2str(model_index),'_',model_type, '.par'];
        md = parameterize(md, ParamFile);
        
        % acquire the sampled var file and substitute the variable
        path2sampled = 'sampled vars/sampled.mat';
        load(path2sampled) % var name is sampled
        
        % get the coordinate data for interpolation in the substitution
        % function
        syn = testbed_data(geometry_path);
        % substitute in
        md = model_param_substitute(md, sampled, syn.x, syn.y);
    
    end
    
    md.transient.ismasstransport=1;
	md.transient.isstressbalance=1;
    md.transient.isgroundingline=1;
	md.transient.ismovingfront=1; % enable levelset method
	md.transient.isthermal=0;
	md.verbose.solution=1;
    md.timestepping.start_time = 0;
    % md.timestepping.time_adapt = 1; % CFL conditon; don't know how to
    % activate
    md.timestepping.time_step = dt;
    if strcmp(model_type, 'spinup')       
        md.timestepping.final_time = nt_spinup*dt;
    else % it is an actual transient run
        md.timestepping.final_time = nt_t*dt;
    end
    
    % Request output
    md.transient.requested_outputs={'default','IceVolume','IceVolumeAboveFloatation'};

    md.cluster=generic('name',oshostname(),'np',n_process);
	md=solve(md,'Transient');

% output
    output = md;
    
%% Saving
% % If spinup run, we save velocity and the model itself
     if strcmp(model_type, 'spinup')
%         vel_grid = griddata(md.mesh.x, md.mesh.y,...
%                                     md.results.TransientSolution(end).Vel,...
%                                     X, Y);
%         vx_grid  = griddata(md.mesh.x, md.mesh.y,...
%                                     md.results.TransientSolution(end).Vx,...
%                                     X, Y);
%         vy_grid  = griddata(md.mesh.x, md.mesh.y,...
%                                     md.results.TransientSolution(end).Vy,...
%                                     X, Y);
%         try % if SSA model, there is no vz         
%             vz_grid  = griddata(md.mesh.x, md.mesh.y,...
%                                         md.results.TransientSolution(end).Vz,...
%                                         X, Y);
%         catch
%             vz_grid  = zeros(size(vy_grid));
%         end
%         
%         V.vel = vel_grid;
%         V.vx  = vx_grid;
%         V.vy  = vy_grid;
%         V.vz  = vz_grid;
%         file_path = ['spinup/spinup_V_',num2str(model_index)];
%         save(file_path,'V')
        
        % save the model, md
        md_file_path = ['spinup_md/spinup_md_', num2str(model_index)];
        save(md_file_path, 'md');
    end
end