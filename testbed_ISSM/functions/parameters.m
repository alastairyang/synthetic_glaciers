function params = parameters(input)
% A callable function which outputs parameters that are of interests
% to this modeling practice
% the input is always a table. If no var is to be changed, input an empty
% table
    %persistent params_persist
        n_layer = 3;
        n_process = 3;
        exponent = 1;
        sim_year_t = 100; % transient run number of iterations
        nt_spinup = 100; % still figuring out if I should specify total iteration number of total sim year...
        sim_year_spinup = 100; % spinup run, number of iterations
        max_stress_grounded = 1000000.0; % 1 MPa
        max_stress_floating = 100000.0; % 100 kPa (0.1*max stress grounded)
        f = 0.6; % Coulomb sliding law coefficient
        % hmin and hmax for the anistropic meshing "bamg"
        hmin = 500;
        hmax = 10000;

        params = table(n_layer, n_process, exponent, sim_year_t, sim_year_spinup, nt_spinup, hmin, hmax, f,...
                       max_stress_grounded, max_stress_floating);

        varnames = input.Properties.VariableNames;
        for i = 1:numel(varnames)
            % update the info
            params.(varnames{i}) = input.(varnames{i});
        end

    %    params_persist = params;
    %else % when the input is not a table
    %    params = params_persist;
    %end
            
end

