function output = my_model_execute_ss(geometry_path, model_index, model_type)
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
    % get a preliminary mesh which we will refine later
    md = bamg(md,'domain', 'Domain.exp', 'hmax', 500); % in meter
    
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
    ParamFile = ['parameters/syn_',model_index,'_ss.par'];
	md = parameterize(md, ParamFile);
    
    %The rheology parameters sit in the material section #md.materials
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

    % set constant smb
    md.smb.mass_balance = InterpFromGridToMesh(syn.x', syn.y, syn.SMB_cons, md.mesh.x, md.mesh.y, 0);
    
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
    vel_grid = griddata(md.mesh.x, md.mesh.y,...
                                    md.results.StressbalanceSolution.Vel,...
                                    syn.X, syn.Y);
    vx_grid  = griddata(md.mesh.x, md.mesh.y,...
                                    md.results.StressbalanceSolution.Vx,...
                                    syn.X, syn.Y);
    vy_grid  = griddata(md.mesh.x, md.mesh.y,...
                                    md.results.StressbalanceSolution.Vy,...
                                    syn.X, syn.Y);
    vz_grid  = griddata(md.mesh.x, md.mesh.y,...
                                    md.results.StressbalanceSolution.Vz,...
                                    syn.X, syn.Y);
    V.vel = vel_grid;
    V.vx  = vx_grid;
    V.vy  = vy_grid;
    V.vz  = vz_grid;
    file_path = ['steady state/ss_V_',model_index];
    save(file_path,'V')
end