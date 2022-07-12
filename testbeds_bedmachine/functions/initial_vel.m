function vel_init = initial_vel(Ug, syn)
%%VEL_INIT generate a synthetic initial velocity distribution using
%%continuity equation.
%  
%   Input:
%       Ug [double]: m/s, grounding line velocity
%       syn[struct]: structure containing geometry data
%
%   Output:
%       vel_init[double array]: velocity field

    %% define parameters
    dx = 150; % m, dx and dy
    dy = dx;
    x_u0_limit = 20000; % m, 10 km before vx drops to zeros (wee need to cap it)

    %% Velocity along the center flowline on the ice shelf
    % Using continuity. A side note that there is analytical solution to ice
    % shelf geometry and velocity field (simultaneous), known as two-part van
    % der Veen exact ice shelf solution (see Ed's note in the textbook folder)
    % the second part of van der Veen exact ice shelf solution is the H(u),
    % which essentially is a continuity relationship. That is all we need as we
    % prescribe the geometry to begin with.
    idx_shelf = find(syn.ocean_mask(:,1) == 0);
    vel_init = zeros(size(syn.s));

    % iterate over each along-flow ice shelf line
    for i = 1:length(idx_shelf)
        % Part A: Ice shelf part
        shelf_line = syn.ocean_mask(idx_shelf(i),:);
        n_shelf = length(shelf_line) - sum(shelf_line); % this is how many elements are marine
        % get the flux
        Hshelf = syn.h(idx_shelf(i), 1:n_shelf);
        Hg = Hshelf(end);
        M0 = Ug * Hg;
        Ushelf = M0./Hshelf;

        % Part B: grounded ice part
        % for the grounded part, we use a hypothetical thickness profile
        % extending from the ice shelf (i.e., linear extrapolation of shelf H).
        % The reason is, there will be a singularity/jump of H(x) and hence
        % U(x) if we adopt the grounded thickness profile, and we don't the
        % singularity to be fed to the initialization v.
        H_slope = (Hg - Hshelf(1))/(length(Hshelf)*dx);
        %Hground_fake = ones(length(shelf_line)-n_shelf, 1);
        x_grounded = syn.X(idx_shelf(i), n_shelf+1:end);
        xg = x_grounded(1) - dx; % x of grounding line
        Hground_fake = Hg + (x_grounded - xg).*H_slope;
        Uground_conti = M0./Hground_fake;
        % another case: velocity linearly decreases to 0 at the grounded
        % ice section
        dudx = Ushelf(end)/(xg - x_u0_limit);
        Uground_linear = Ushelf(end) + (x_grounded - xg).*dudx;
        Uground_linear(Uground_linear < 0) = 0;
        Uground = zeros(1, size(syn.s, 2)-n_shelf);
        Uground(Uground_linear<Uground_conti) = Uground_linear(Uground_linear<Uground_conti);
        Uground(Uground_linear>=Uground_conti) = Uground_conti(Uground_linear>=Uground_conti);

        % Assign to the vel array
        vel_row = [Ushelf, Uground];
        vel_init(idx_shelf(i),:) = vel_row;
        % although the thickness is a linear extrapolation of shelf thickness,
        % thanks to the continuity equation, the velocity field still
        % diminished inversely wrt distance from the calving front
    end
    % smooth the filed
    vel_init = imgaussfilt(vel_init,'FilterSize',3);
end

    
