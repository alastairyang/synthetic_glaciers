function output = t_my_model_execute_ss(geometry_path, model_index, model_type)
%%  This function excute the model without calling the parameter files
%%  Model
%   As the parameter files will be called in by parameters
    md = model;
    % parse out params
    % calling parameters(1) is to input something other than a table
    % so nothing is changed (changes should be made in the master file)
    % and hence this function is just to output
    
    % create an empty table for storing
    tb = table();
    % the parameters() function is simply to store hyperparamters of
    % this model; they are not physical parameters
    params = parameters(tb);
    n_layer = params.n_layer;
    n_process = params.n_process;
    exponent = params.exponent;
	
    % we use triangle method to create mesh
    % Thi is modeled after JI example
    % triangle(model, domain_file, average_element_size_meter)
    syn = testbed_data(geometry_path);
    [Xq,Yq,~] = meshgrid_downsample(syn.X, syn.Y, ones(size(syn.X))); % use ones as placeholder
    
    % generate Domain.exp file
    meshgrid2outline(Xq,Yq);
    % get a mesh
    md = bamg(md,'domain', 'Domain.exp', 'hmax', 1500); % in meter
    
    %md=triangle(md, 'Domain.exp');
    % plotmodel(md, 'data','mesh')
    
    % refine the mesh WRT a field data
    % first, check what the vel looks like
    % plotmodel(md,'data',vel,'edgecolor','w')
    % md = bamg(md, 'domain','Domain.exp','field', vel, 'err', 0.05);
    
%Masks #2
    groundedice = double(InterpFromGridToMesh(syn.x',syn.y,syn.ocean_mask,md.mesh.x,md.mesh.y,0));
	groundedice(groundedice <= 0) = -1; % when 0, make it negative, hence no ice
    icepresence = double(InterpFromGridToMesh(syn.x',syn.y,syn.ice_mask  ,md.mesh.x,md.mesh.y,0));

    %fill in the md.mask structure
	md.mask.ocean_levelset=groundedice; %ice is grounded for mask equal one
	md.mask.ice_levelset=icepresence;%ice is present when negatvie, and it is present everywhere

	%ploting
	plotmodel(md,'data', 'mesh', 'data',md.mask.ocean_levelset,'title','grounded/floating','data',md.mask.ice_levelset,'title','ice/no-ice')
	
    %Parameterization #3
%     ParamFile = ['parameters/syn_',model_index,'_ss.par'];
% 	md = parameterize(md, ParamFile);

    %Geometry
    disp('   Constructing Geometry');

    %Define the geometry of the simulation #md.geometry
    % S: surface elevation
    % B: base elevation
    % H: ice thickness
    % X, Y: Nx1 coordinates
    % X_g, Y_g: meshgrid
    disp('   Loading geometry data');
    [geometry, ~, ~, ~] = query_data(model_index,model_type);
    geometry_path = geometry{1};
    syn = testbed_data(geometry_path);
    md.geometry.surface = InterpFromGridToMesh(syn.x', syn.y, syn.s, md.mesh.x, md.mesh.y, 0);
    md.geometry.bed  = InterpFromGridToMesh(syn.x', syn.y, syn.bed,  md.mesh.x, md.mesh.y, 0);
    md.geometry.base = md.geometry.bed;
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

    disp('   Defining friction parameters');
    fric_coef_mesh = InterpFromGridToMesh(syn.x', syn.y, syn.fric_coef,...
                                          md.mesh.x, md.mesh.y, 0);
    md.friction.coefficient=fric_coef_mesh;
    %one friciton exponent (p,q) per element
    md.friction.p=ones(md.mesh.numberofelements,1);
    md.friction.q=ones(md.mesh.numberofelements,1);

    %no friction applied on floating ice
    pos=find(md.mask.ocean_levelset<0);
    md.friction.coefficient(pos)=0;
    md.groundingline.migration='SubelementMigration';

    disp('   Construct ice rheological properties');

    % %The rheology parameters sit in the material section #md.materials
    % %B has one value per vertex
    % md.materials.rheology_B=6.8067e7*ones(md.mesh.numberofvertices,1);
    %n has one value per element
    %->
    md.materials.rheology_n=3*ones(md.mesh.numberofelements,1);

    disp('   Set boundary conditions');
    md=SetMarineIceSheetBC(md);

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
    % also fix the thickness at the calving front
    % well, basically fixing all boundary thickness
    boundary_nodes = find(md.mesh.vertexonboundary);
    md.masstransport.spcthickness(boundary_nodes) = md.geometry.thickness(boundary_nodes);
    
    %B has one value per vertex
    
    md.materials.rheology_B = double(InterpFromGridToMesh(syn.x',syn.y,syn.rheoB.rheoB_unif,md.mesh.x,md.mesh.y,0));
    
    %Extrusion #?
	% only 5 layers exponent 1
	md = extrude(md, n_layer, exponent);
    
    %Forcing: needs to be forced after extrusion for model consistency
%     % Ice shelf basal forcing
%     shelf_melt_mesh = InterpFromGridToMesh(x', y, shelf_melt.constant_melt, md.mesh.x, md.mesh.y, 0);
%     md.basalforcings.floatingice_melting_rate = shelf_melt_mesh;
% 
%     % set the surface mass balance
%     md.smb = SMBgradients();
%     href_mesh   = InterpFromGridToMesh(x', y, SMB_grad.href, md.mesh.x, md.mesh.y, 0);
%     smbref_mesh = InterpFromGridToMesh(x', y, SMB_grad.smbref, md.mesh.x, md.mesh.y, 0);
%     md.smb.href   = href_mesh;
%     md.smb.smbref = smbref_mesh;
%     % b_pos and b_neg also need to be fields
%     md.smb.b_pos = repmat(SMB_grad.b_pos, md.mesh.numberofvertices,1); 
%     md.smb.b_neg = repmat(SMB_grad.b_neg, md.mesh.numberofvertices,1);
    
    % set name
    md.miscellaneous.name = ['ss_',model_index];
    
%Set the flow computing method #5
	md = setflowequation(md, 'HO', 'all');
    % before running the model, take a look at all boundary conditions
    plotmodel(md,'data','BC')

%Solving #6
	md.cluster=generic('name',oshostname(),'np',n_process);
	md.verbose=verbose('convergence',true);

	md=solve(md,'Stressbalance');

% output
    output = md;
    
% plot
    plotmodel(md, 'data', md.results.StressbalanceSolution.Vel)
    
% save the velocity field data
    % use the InterpFromMeshToGrid tool
    vel_grid = InterpFromMeshToGrid(md.mesh.elements, md.mesh.x, md.mesh.y,...
                                    md.results.StressbalanceSolution.Vel,...
                                    syn.x, syn.y);
    vx_grid  = InterpFromMeshToGrid(md.mesh.elements, md.mesh.x, md.mesh.y,...
                                    md.results.StressbalanceSolution.Vx,...
                                    syn.x, syn.y);
    vy_grid  = InterpFromMeshToGrid(md.mesh.elements, md.mesh.x, md.mesh.y,...
                                    md.results.StressbalanceSolution.Vy,...
                                    syn.x, syn.y);
    vz_grid  = InterpFromMeshToGrid(md.mesh.elements, md.mesh.x, md.mesh.y,...
                                    md.results.StressbalanceSolution.Vz,...
                                    syn.x, syn.y);
    V.vel = vel_grid;
    V.vx  = vx_grid;
    V.vy  = vy_grid;
    V.vz  = vz_grid;
    file_path = ['steady state/ss_V_',model_index];
    save(file_path,'V')
end