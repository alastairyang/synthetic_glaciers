
%Geometry
disp('   Acquire data');
[geometry, ~] = query_data(3,'spinup');
geometry_path = geometry{1}; % this is actually only for importing fric_coef later

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
[x, y, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, fric_coef, ~, ~] = testbed_data(geometry_path);
fric_coef_mesh = InterpFromGridToMesh(x', y, fric_coef, md.mesh.x, md.mesh.y, 0);
md.friction.coefficient = fric_coef_mesh;
pos = find(md.mask.ocean_levelset<0);
md.friction.coefficient(pos) = 0;

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

disp('    Initiate levelset')
md.levelset.spclevelset = md.results.TransientSolution(end).MaskIceLevelset;
icefront_i = find(vertexonicefront == 1);
md.levelset.spclevelset(icefront_i) = 0; % ice front needs to be zero


disp('    Done importing the spinup model info; clear the field')
% clear the solution field
md.results = rmfield(md.results, 'TransientSolution');