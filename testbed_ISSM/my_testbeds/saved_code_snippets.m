%% Scenario: how to better mesh the spin-up runs
% the idea was that
%    1. since stress balance solution always get it wrong, generating a
%    legit velocity field is untenable that way.
%    2. Instead we can run a transient for a few years and then use the
%    velocity fields from that to re-mesh the spinup.
%    3. We have a working fine mesh without field adaptation, but it is
%    fine everywhere and that lengthens the computational time.
% Why is this abandoned:
%    1. in bamg, 'err' controls significantly the mesh gradation. With the
%       default error, which we don't know the value, the mesh is fine
%       everywhere. Kinda meaningless.
        md_spinup_prelim = load(['spinup prelim/',model_index,'.mat']);
        md_spinup_prelim = md_spinup_prelim.md;
        prelim_vx = md_spinup_prelim.results.TransientSolution(end).Vx;
        prelim_vy = md_spinup_prelim.results.TransientSolution(end).Vy;
        prelim_vx_grid = InterpFromMeshToGrid(md_spinup_prelim.mesh.elements,...
                                              md_spinup_prelim.mesh.x, md_spinup_prelim.mesh.y, prelim_vx,...
                                              syn.x, syn.y, mean(prelim_vx,'all'));
        prelim_vy_grid = InterpFromMeshToGrid(md_spinup_prelim.mesh.elements,...
                                              md_spinup_prelim.mesh.x, md_spinup_prelim.mesh.y, prelim_vy,...
                                              syn.x, syn.y, mean(prelim_vy,'all'));                        
        prelim_vx_mesh = InterpFromGridToMesh(syn.x', syn.y, prelim_vx_grid, md.mesh.x, md.mesh.y, mean(prelim_vx_grid,'all'));                                   
        prelim_vy_mesh = InterpFromGridToMesh(syn.x', syn.y, prelim_vx_grid, md.mesh.x, md.mesh.y, mean(prelim_vy_grid,'all'));                                   
        
        prelim_vel = sqrt(prelim_vx_mesh.^2 + prelim_vy_mesh.^2);
        prelim_vel_trans = 500*(normalize(prelim_vel) - min(normalize(prelim_vel)));
        md = bamg(md,'hmin',500,'hmax',200000,'field',prelim_vel_trans, 'err',20);