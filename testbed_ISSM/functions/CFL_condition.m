function dt = CFL_condition(vel, dx_min, dy_min)
%CFL_CONDITION CFL condition. We use the two components of velocity field
%from steady-state simulation results to calculate the delta t (time-step)
%from CFL condtion, assuming a Courant number.
%
%   Input:
%       vel_fields: structure, storing the two velocity components
%       dx_min: double, the min dx from the "bamg" meshgrid
%       dy_min: double, the min dy from the "bamg" meshgrid

    % parameter: Courant number
    Cmax = 0.6;
    
    % velocities
    % since we are using SSA, only x and y components
    vx_max = max(vel, [], 'all');
    vy_max = 0;
    % CFL condition on two dimensions
    dt = Cmax/(vx_max/dx_min + vy_max/dy_min);
    
end

