function output = my_model_execute_t(geometry_path, velocity_path, model_path, model_index, model_type, forcing)
%% get global variable(s)
global to_disk
%%  Model 
    tb = table();
    % the parameters() function is simply to store hyperparamters of
    % this model; they are not physical parameters
    params = parameters(tb);
    n_layer = params.n_layer;
    n_process = params.n_process;
    exponent = params.exponent;
    sim_year_t = params.sim_year_t;
    sim_year_spinup = params.sim_year_spinup;
    nt_spinup = params.nt_spinup;
    max_stress_grounded = params.max_stress_grounded;
    max_stress_floating = params.max_stress_floating;
    hmin = params.hmin;
    hmax = params.hmax;
    f = params.f; % the f in Coulomb friction law

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
        [Xq,Yq,~] = meshgrid_downsample(syn.X, syn.Y, ones(size(syn.X)), model_type); % use ones as placeholder

        % generate Domain.exp file
        meshgrid2outline(Xq,Yq);
        % get a preliminary mesh which we will refine later
        md = bamg(md,'domain', 'Domain.exp', 'hmax', 200);

        % Velocity field
        % load(velocity_path); % loaded as 'V'
        vel_mesh = InterpFromGridToMesh(syn.x',syn.y, syn.vel_init_x,...
                                        md.mesh.x,md.mesh.y, mean(syn.vel_init_x,'all'));
        vel_mesh_norm = normalize(vel_mesh) + min(normalize(vel_mesh), [], 'all');
        md = bamg(md,'hmin',hmin,'hmax',hmax,'field',vel_mesh_norm,'err',5);
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
        plotmodel(md,'data','mesh', 'title','mesh',...
                     'data',md.mask.ocean_levelset,'title','grounded/floating')

    %%  Parameterization #3
%         ParamFile = ['parameters/syn_',model_index,'_',model_type, '.par'];
%         md = parameterize(md, ParamFile);
        disp('   Loading geometry data');
        [geometry, ~] = query_data(model_index,'spinup');
        geometry_path = geometry{1};
        % velocity_path = velocity{1};
        syn = testbed_data(geometry_path);
        md.geometry.surface = InterpFromGridToMesh(syn.x', syn.y, syn.s, md.mesh.x, md.mesh.y, 0);
        md.geometry.bed  = InterpFromGridToMesh(syn.x', syn.y, syn.bed,  md.mesh.x, md.mesh.y, 0);
        md.geometry.base = md.geometry.bed; % we will change this later, the part where ice is floating
        md.geometry.thickness = md.geometry.surface - md.geometry.base;

        %Get the node numbers of floating nodes
        pos=find(md.mask.ocean_levelset<0); 

        %apply a flotation criterion on the precedingly defined nodes and
        %redefine base and thickness accordingly
        %ensure hydrostatic equilibrium on ice shelf: 
        di=md.materials.rho_ice/md.materials.rho_water;
        md.geometry.thickness(pos)=1/(1-di)*md.geometry.surface(pos);
        md.geometry.base(pos)=md.geometry.surface(pos)-md.geometry.thickness(pos);

        disp('   Constructing thickness related variables');
        pos0=find(md.geometry.thickness<=1);
        md.geometry.thickness(pos0)=1;
        md.geometry.surface=md.geometry.thickness+md.geometry.base;
        md.geometry.hydrostatic_ratio=ones(md.mesh.numberofvertices,1);


        disp('   Defining friction parameters');
        fric_coef_mesh = InterpFromGridToMesh(syn.x', syn.y, syn.fric_coef,...
                                              md.mesh.x, md.mesh.y, 0);
        md.friction = frictiontsai();
        md.friction.C=fric_coef_mesh;
        % the m in Weertman's/Tsai's laws is the inverse of p in Paterson's formulation
        % md.friction.m = syn.slidingP*ones(md.mesh.numberofelements,1);
        md.friction.m = (1/syn.slidingP)*ones(md.mesh.numberofelements, 1);
        md.friction.f = f*ones(md.mesh.numberofelements, 1);
        
        %no friction applied on floating ice
        pos=find(md.mask.ocean_levelset<0);
        md.friction.C(pos)=0.0000000001;
        md.groundingline.migration='SubelementMigration';

        disp('   Construct ice rheological properties');

        % %The rheology parameters sit in the material section #md.materials
        % %B has one value per vertex
        % md.materials.rheology_B=6.8067e7*ones(md.mesh.numberofvertices,1);
        %n has one value per element
        %->
        md.materials.rheology_n=3*ones(md.mesh.numberofelements,1);

        disp('   Set boundary conditions');

        %Set the default boundary conditions for an ice-sheet 
        % #help SetIceSheetBC
        md=SetMarineIceSheetBC(md);

        % Initializing: pressure, velocity field
%         load(velocity_path) % loaded as 'V'
%         md.initialization.pressure=md.materials.rho_ice*md.constants.g*md.geometry.thickness;
%         vx_mesh = InterpFromGridToMesh(syn.x', syn.y, V.vx, md.mesh.x, md.mesh.y, 0);
%         vy_mesh = InterpFromGridToMesh(syn.x', syn.y, V.vy, md.mesh.x, md.mesh.y, 0);
%         vz_mesh = InterpFromGridToMesh(syn.x', syn.y, V.vz, md.mesh.x, md.mesh.y, 0);
        vx_mesh = InterpFromGridToMesh(syn.x', syn.y, syn.vel_init_x, md.mesh.x, md.mesh.y, 0);
        vy_mesh = zeros(size(vx_mesh));
        vz_mesh = zeros(size(vx_mesh));
        md.initialization.vx=vx_mesh;
        md.initialization.vy=vy_mesh;
        md.initialization.vz=vz_mesh;
        md.initialization.vel=sqrt(md.initialization.vx.^2+...
                                   md.initialization.vy.^2+...
                                   md.initialization.vz.^2);

        disp('   Fixing thickness at inflow boundary for spin-up run');
        % Set Dirichlet B.C of thickness at inflow boundary
        % first get vertex on ice front (borrowed from SetMarineIceSheet)
        % then subtract
        pos = find(sum(md.mask.ocean_levelset(md.mesh.elements)<0.,2) >0.);
        vertexonfloatingice = zeros(md.mesh.numberofvertices,1);
        vertexonfloatingice(md.mesh.elements(pos,:)) = 1.;
        vertexonicefront = double(md.mesh.vertexonboundary & vertexonfloatingice);
        vertexnoticefront = md.mesh.vertexonboundary - vertexonicefront;
        noticefront_i  = find(vertexnoticefront == 1);
        md.masstransport.spcthickness(noticefront_i) = md.geometry.thickness(noticefront_i);
        
        disp('    Initiate levelset')
        icepresence=double(InterpFromGridToMesh(syn.x',syn.y,syn.ice_mask, md.mesh.x,md.mesh.y,0));
        md.levelset.spclevelset = icepresence;
        icefront_i = find(vertexonicefront == 1);
        md.levelset.spclevelset(icefront_i) = 0; % ice front needs to be zero

    %%  Processes
        
        % rheology - B: data{1} should be uniform B, no shear margin
        % weakening. See it is using the first cell (no weakening)
        md.materials.rheology_B = double(InterpFromGridToMesh(syn.x',syn.y,syn.transient_rheoB.data{1},...
                                                              md.mesh.x,md.mesh.y, mean(syn.transient_rheoB.data{1},'all')));

        % friction coef
        % no perturbation (slippery patch)
        md.friction.C = double(InterpFromGridToMesh(syn.x',syn.y,syn.transient_fric_coef.data{1},...
                                                              md.mesh.x,md.mesh.y, mean(syn.transient_fric_coef.data{1},'all')));
        % no friction applied on floating ice; we did this before in .par
        % but a little redundancy does not hurt
        pos=find(md.mask.ocean_levelset<0);
        md.friction.C(pos)=0.0000000001;
        
        % forcing: needs to be forced after extrusion for model consistency
        % Ice shelf basal melt: constant melt rate
        % Constant melt rate forcing
        shelf_melt_mesh = InterpFromGridToMesh(syn.x', syn.y, syn.shelf_melt.transient_melt{1}, md.mesh.x, md.mesh.y, 0);
        md.basalforcings.floatingice_melting_rate = shelf_melt_mesh;
        
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
        load(model_path) % the variable is md
        
        %% Parameterize
        % This .par substitutes much of the final state of the spin-up run
        % to the initial conditions of this transient run
%         ParamFile = ['parameters/syn_',model_index,'_',model_type, '.par'];
%         md = parameterize(md, ParamFile);
        disp('   Acquire data');
        [geometry, ~] = query_data(model_index,'spinup');
        geometry_path = geometry{1};
        disp('   Substitute the geometry')
        surface  = md.results.TransientSolution(end).Surface;
        base     = md.results.TransientSolution(end).Base;
        thickness = surface - base;

        md.geometry.surface   = surface;
        md.geometry.base      = base;
        md.geometry.thickness = thickness;
        % no changes to md.geometry.bed
        % substitute the mask and redefine 0 friction
        md.mask.ocean_levelset = md.results.TransientSolution(end).MaskOceanLevelset;

        disp('    Min ice thickness is 1 meter!')
        pos0=find(md.geometry.thickness<=1);
        md.geometry.thickness(pos0)=1;
        md.geometry.surface=md.geometry.thickness+md.geometry.base;

        disp('    Substitute friciton coefficent')
        syn = testbed_data(geometry_path);
        fric_coef_mesh = InterpFromGridToMesh(syn.x', syn.y, syn.fric_coef, md.mesh.x, md.mesh.y, 0);
        md.friction.C = fric_coef_mesh;
        pos = find(md.mask.ocean_levelset<0);
        md.friction.C(pos) = 0;

        disp('    Substitute the initial velocity distribution')
        % get the velocity from results and put it into initialization
        Vel = md.results.TransientSolution(end).Vel;
        Vx  = md.results.TransientSolution(end).Vx;
        Vy  = md.results.TransientSolution(end).Vy;
        try % if SSA, there is no vz, then just 0 all
            Vz  = md.results.TransientSolution(end).Vz;
        catch
            Vz  = zeros(size(Vy));
        end
        md.initialization.vx  = Vx;
        md.initialization.vy  = Vy;
        md.initialization.vz  = Vz;
        md.initialization.vel = Vel;

        disp('   Constructing boundary conditions')
        % clear the spcthickness condition
        md.masstransport.spcthickness = NaN*ones(size(md.geometry.thickness));

        % Vertex on ice front (borrowed from SetMarineIceSheet)
        pos = find(sum(md.mask.ocean_levelset(md.mesh.elements)<0.,2) >0.);
        vertexonfloatingice = zeros(md.mesh.numberofvertices,1);
        vertexonfloatingice(md.mesh.elements(pos,:)) = 1.;
        vertexonicefront = double(md.mesh.vertexonboundary & vertexonfloatingice);
        vertexnoticefront = md.mesh.vertexonboundary - vertexonicefront;
        % Get the thickness and fix the thickness there
        fixh_i = find(vertexnoticefront == 1);
        md.masstransport.spcthickness(fixh_i) = md.geometry.thickness(fixh_i);

        disp('    Initialize levelset')
        md.levelset.spclevelset = md.results.TransientSolution(end).MaskIceLevelset;
        icefront_i = find(vertexonicefront == 1);
        md.levelset.spclevelset(icefront_i) = 0; % ice front needs to be zero

        %% Load forcing
        disp('    Done importing the spinup model info; clear the field')
        % clear the solution field
        md.results = rmfield(md.results, 'TransientSolution');
        
        % get data
        syn = testbed_data(geometry_path);
        
        if ~isempty(forcing)
            if strcmp('all',forcing)
                disp('Forcing: all forcing, ON')
                pause(1)
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
                        md.friction.C = temp;
                    else
                        md.friction.C = [md.friction.C, temp];
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
                
            elseif strcmp('fric',forcing)
                disp('Forcing: fric.coef perturbation, ON')
                pause(1)
                % time-dependent frictional coefficient
                N_years = numel(syn.transient_fric_coef.years);
                for i = 1:N_years
                    fric_coef_mesh = InterpFromGridToMesh(syn.x', syn.y, syn.transient_fric_coef.data{i}, md.mesh.x, md.mesh.y, 0);
                    temp = [fric_coef_mesh; syn.transient_fric_coef.years{i}]; % vertical concat year
                    if i == 1
                        md.friction.C = temp;
                    else
                        md.friction.C = [md.friction.C, temp];
                    end
                end
                
            elseif strcmp('rheoB',forcing)
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
                
            elseif strcmp('meltrates',forcing)
                disp('Forcing: submarine melting, ON')
                pause(1)
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
            else
                disp('No forcing? Strings do not match')
                pause(2)
                return
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
        ParamFile = ['parameters/syn_',model_index,'_',model_type, '.par'];
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
    md.transient.ismovingfront=1;
    md.transient.isgroundingline=1;
	md.transient.isthermal=0;
	md.verbose.solution=1;
    md.timestepping.start_time = 0;
    % timestep, smaller one between 0.1 yr and dt from CFL condition
    md.timestepping.time_step = min(CFL_condition(syn.vel_init_x, hmin, hmin), 0.1);
    disp(['Time step is ', num2str(md.timestepping.time_step)])
    if strcmp(model_type, 'spinup')       
        md.timestepping.final_time =  md.timestepping.time_step*nt_spinup;
    else % it is an actual transient run
        md.timestepping.final_time = sim_year_t;
    end
    
    % Request output
    md.transient.requested_outputs={'default','IceVolume','IceVolumeAboveFloatation'};

    md.cluster=generic('name',oshostname(),'np',n_process);
    md.miscellaneous.name = 'transient';
	md=solve(md,'Transient');

% output
    output = md;
    
%% Saving
% % If spinup run, we save the model itself
     if strcmp(model_type, 'spinup')

        % save the model, md
        if to_disk % global variable.
            md_file_path = ['/Volumes/Donglai_SSD/spinup_md/spinup_md_', model_index];
        else
            md_file_path = ['spinup_md/spinup_md_', model_index];
        end
        save(md_file_path, 'md');
        % and save some plots to monitor if relaxation is completed
        export_graph(md, model_index, model_type)
        
     elseif strcmp(model_type, 't')
        foldername = ['syn_', model_index];
        md_file_path = ['results/', foldername, '/t_md_', model_index,'_', forcing];
        save(md_file_path, 'md');
    end
end