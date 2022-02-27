function [geo_filepath, v_filepath, md_filepath, sens_filepath] = query_data(model_index,model_type)
%QUERY_DATA Retrieve the path to geometry or geometry/velocity data file
%depending on the model type
%
%   Input:
%       model_index[string]: model index
%       model_type[string] : either steady-state or transient
%
%   Output:
%       data_paths[struc]  : a structure storing the paths to data files
    current_dir = pwd;
    % get the path of the gemoetry file given the model index
    geo_filename = ['syn_', model_index,'.mat'];
    geo_filepath = get_absolute_path('git_research',geo_filename);
    v_filepath  = [];
    md_filepath = [];
    sens_filepath = []; % if type is 'ss', these and only these lines are done
    
    if strcmp(model_type, 'spinup') % spinup
        % get the velocity file
        v_filename = ['ss_V_', model_index,'.mat'];
        v_filepath = get_absolute_path('git_research',v_filename);
        md_filepath = [];
        sens_filepath = [];
    end
    
    if strcmp(model_type, 't') % actual transient forced with time-dependent climate
        % get the velocity file from spinup with fixed thickness
        v_filename = ['spinup_V_', model_index,'.mat'];
        v_filepath = get_absolute_path('git_research',v_filename);
        md_filename = ['spinup_md_', model_index,'.mat'];
        md_filepath = get_absolute_path('git_research',md_filename);
        sens_filepath = [];
    end

    if strcmp(model_type, 't_sensitive') % Runs with perturbed transient forcings
        % get the velocity file from spinup with fixed thickness
        v_filename = ['spinup_V_', model_index,'.mat'];
        v_filepath = get_absolute_path('git_research',v_filename);
        md_filename = ['spinup_md_', model_index,'.mat'];
        md_filepath = get_absolute_path('git_research',md_filename);
        sens_filename = ['sens_', model_index, '.mat'];
        sens_filepath = get_absolute_path('git_research', sens_filename);
    end
    
    cd(current_dir)
end

