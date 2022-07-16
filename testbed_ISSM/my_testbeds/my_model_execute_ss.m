function my_model_execute_ss(geometry_path, model_index, model_type)
%%  This function excute the model without calling the parameter files
%%
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
        md = bamg(md,'domain', 'Domain.exp', 'hmax', syn.attrs_table.fjord_width*0.3, 'hmin',100, 'gradation',1.01);
        
        plotmodel(md, 'data','mesh')

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
        disp('   Make min thickness 1 meter');
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
        %md=SetMarineIceSheetBC(md);

        % Initializing: pressure, velocity field
%         load(velocity_path) % loaded as 'V'
%         vx_mesh = InterpFromGridToMesh(syn.x', syn.y, V.vx, md.mesh.x, md.mesh.y, 0);
%         vy_mesh = InterpFromGridToMesh(syn.x', syn.y, V.vy, md.mesh.x, md.mesh.y, 0);
%         vz_mesh = InterpFromGridToMesh(syn.x', syn.y, V.vz, md.mesh.x, md.mesh.y, 0);
        disp('    Initialization');
        md.initialization.pressure=md.materials.rho_ice*md.constants.g*md.geometry.thickness;
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
        
        % Set boundary condition
        md=SetMarineIceSheetBC(md);
        
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
        
        md.transient.ismasstransport=1;
        md.transient.isstressbalance=1;
        md.transient.ismovingfront=1;
        md.transient.isgroundingline=1;
        md.transient.isthermal=0;
        md.verbose.solution=1;
        md.timestepping=timesteppingadaptive();
        md.timestepping.start_time = 0;
        md.timestepping.final_time =  5;


        % Request output
        md.transient.requested_outputs={'default','IceVolume','IceVolumeAboveFloatation'};

        md.cluster=generic('name',oshostname(),'np',n_process);
        md.miscellaneous.name = model_type;

        %md = solve(md,'Stressbalance');
        md = solve(md,'Transient');
        % output
        
        file_path = ['spinup prelim/',model_index,'.mat'];
        save(file_path,'md')
end