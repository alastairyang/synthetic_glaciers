function output_model = run_models(model_index, model_type, forcing)
%RUN_MODELS Control what model to run
%
%   Input:
%       model_index[string]: model index
%       model_type[string] : either steady-state or transient
%       forcing[string]    : a specific variable/forcing that is activated
%  
%   Output:
%       output_model[struc]: a structure storing the ISSM model outputs
    
    model_name   = ['model_', model_index];

    if strcmp(model_type, 'spinup') % spin up
        % get the velocity and geometry file
        [geometry, ~, ~, ~] = query_data(model_index, model_type);
        % input to model
        md_spinup = my_model_execute_t(geometry{1}, [], [], model_index, model_type,[]);
        output_model.([model_name,'_spinup']) = md_spinup;
    end
    
    if strcmp(model_type, 't') % actual transient run forced with time-dependent climate
        % get the velocity and geometry file
        [geometry, ~, model, ~] = query_data(model_index, model_type);
        % input to model
        md_t = my_model_execute_t(geometry{1}, [], model{1}, model_index, model_type, forcing);
        output_model.([model_name,'_t']) = md_t;  
    end
    
    if strcmp(model_type, 't_sensitive')
        % on top of the 't' run, we need to get the sensitivity data
        % get the velocity and geometry file
        [geometry, velocity, model, sens_data] = query_data(model_index, model_type,forcing);
        sens_data = load(sens_data{1});
        dataname = ['sens_',num2str(model_index)];
        N_test = sens_data.(dataname).N_test;
        N_vars  = sens_data.(dataname).N_vars;
        for test_index = 1:N_test*N_vars
            % Sample a combination of perturbed forcings; this should
            % generate a file in 'sampled vars' folder. No output for here.
            sens_data_sampler(model_index, test_index, sens_data)
            % each sensitivity test iteration is performed here
            md_t = my_model_execute_t(geometry{1}, velocity{1}, model{1}, model_index, model_type,[]);
            output_model.([model_name,'_sens_test',num2str(test_index)]) = md_t; 
        end
    end
    
    if strcmp(model_type, 'ss') % transient
        [geometry, ~, ~, ~] = query_data(model_index, model_type);
        % input to model
        my_model_execute_ss(geometry{1}, model_index, model_type);
        % output_model.([model_name,'_ss']) = md_ss;
    end
    
end

