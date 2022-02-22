function output = model_execute_t()
%%  This function excute the model without calling the parameter files
%   As the parameter files will be called in by parameters
    md = model;
    % parse out params
    % calling parameters(1) is to input something other than a table
    % so nothing is changed (changes should be made in the master file)
    % and hence this function is just to output
    tb = table();
    % the parameters() function is simply to store hyperparamters of
    % this model; they are not physical parameters
    params = parameters(tb);
    n_layer = params.n_layer;
    n_process = params.n_process;
    exponent = params.exponent;
    dt = params.dt;
    nt = params.nt;
	
    % we use triangle method to create mesh
    % Thi is modeled after JI example
    % triangle(model, domain_file, average_element_size_meter)
    [x, y, X, Y,~,~,~, VX, VY] = testbed_data();
    [Xq,Yq,~] = meshgrid_downsample(X, Y, ones(size(X))); % use ones as placeholder
    
    % generate Domain.exp file
    meshgrid2outline(Xq,Yq);
    % get a preliminary mesh which we will refine later
    md = bamg(md,'domain', 'Domain.exp', 'hmax', 3000);
    
    % Velocity field
    VX_m = InterpFromGridToMesh(x,y,VX,md.mesh.x,md.mesh.y,0);
	VY_m = InterpFromGridToMesh(x,y,VY,md.mesh.x,md.mesh.y,0);
    % Speed
	vel = sqrt(VX_m.^2+VY_m.^2);
    md = bamg(md,'hmin',800,'hmax',15000,'field',vel,'err',5);
    %md=triangle(md, 'Domain.exp');
    % plotmodel(md, 'data','mesh')
    
    % refine the mesh WRT a field data
    % first, check what the vel looks like
    % plotmodel(md,'data',vel,'edgecolor','w')
    % md = bamg(md, 'domain','Domain.exp','field', vel, 'err', 0.05);
    
%Masks #2
	% all MISMIP nodes are grounded
	md = setmask(md, '', '');
    
%Parameterization #3
    ParamFile = 'testbed_t.par';
	md = parameterize(md, ParamFile);

%Extrusion #4
	% only 5 layers exponent 1
	md = extrude(md, n_layer, exponent);

%Set the flow computing method #5
	md = setflowequation(md, 'HO', 'all');

%Set Boundary Conditions #6
% 	md.stressbalance.spcvx(:) = NaN*ones(md.mesh.numberofvertices, 1);
%     md.stressbalance.spcvy(:) = NaN*ones(md.mesh.numberofvertices, 1);
%     md.stressbalance.spcvz(:) = NaN*ones(md.mesh.numberofvertices, 1);
% 	% extract the nodenumbers at the base #md.mesh.vertexonbase
% 	basalnodes = find(md.mesh.vertexonbase);
%   % set the sliding to zero on the bed (Vx and Vy)
% 	md.stressbalance.spcvx(basalnodes)=0.0;
% 	md.stressbalance.spcvy(basalnodes)=0.0;
% 	% periodic boundaries have to be fixed on the sides
% 	maxX = find(md.mesh.x == max(md.mesh.x));
% 	minX = find(md.mesh.x == min(md.mesh.x));
% 	% for y, max X and minX should be excluded
% 	maxY = find(md.mesh.y == max(md.mesh.y) & md.mesh.x ~= max(md.mesh.x) & md.mesh.x ~= min(md.mesh.x));
% 	minY = find(md.mesh.y == min(md.mesh.y) & md.mesh.x ~= max(md.mesh.x) & md.mesh.x ~= min(md.mesh.x));
% 	md.stressbalance.vertex_pairing = [minX, maxX; minY, maxY];

%Solving #8
	md.cluster=generic('name',oshostname(),'np',n_process);
	md.verbose=verbose('convergence',true);

    md.transient.ismasstransport=1;
	md.transient.isstressbalance=1;
	md.transient.ismovingfront=0;
	md.transient.isthermal=0;
	md.verbose.solution=1;
    md.timestepping.time_step = dt;
	md.timestepping.final_time = dt*nt;
    
	md=solve(md,'Transient');

% output
    output = md;
    
% plot
    plotmodel(md,'data', md.results.TransientSolution(1).Vel,...
                 'data', md.results.TransientSolution(nt).Vel)
end