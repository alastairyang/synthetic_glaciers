function params = parameters(input)
% A callable function which outputs parameters that are of interests
% to this modeling practice
% the input is always a table. If no var is to be changed, input an empty
% table
    %persistent params_persist
    
    %if isequal(class(input), 'table')
        n_layer = 3;
        n_process = 3;
        exponent = 1;
        dt = 0.1; % 
        nt_t = 1000; % transient run number of iterations
        nt_spinup = 1000; % spinup run, number of iterations
        max_stress_grounded = 1000000.0; % 1 MPa
        max_stress_floating = 100000.0; % 100 kPa (0.1*max stress grounded)

        params = table(n_layer, n_process, exponent, dt, nt_t, nt_spinup, max_stress_grounded, max_stress_floating);

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

