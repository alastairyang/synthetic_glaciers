%% Spin-up all the geometries
clear; clc;
cd '/Users/donglaiyang/Library/Mobile Documents/com~apple~CloudDocs/Documents/Documents_Alastair’s_MacBook_Pro/Buffalo/Research/git_research/testbed_ISSM/my_testbeds'
all_filenames = struct2table(dir('../../testbeds_bedmachine/Synthetic glaciers/newlabels/*.mat'));
name_split = split(all_filenames.name, '_');
model_labels = split(name_split(:,2), '.');
model_labels = model_labels(:,1);
geom_paths = all_filenames.folder;

% we only do spinup runs for this script
model_type = 'spinup';
global to_disk
to_disk = true;

%% Run the models
N_model = length(model_labels);
for i = 1:N_model
    path = geom_paths(i);
    model_label = model_labels(i);
    model_label = model_label{:};
    progress_count = ['(', num2str(i), '/', num2str(N_model), ')'];
    disp(['Model ', model_label, ' is starting, ', progress_count])
    
    % first, run the steady-state model
    run_models(model_label, 'ss',       []);
    % then, run the spinup
    run_models(model_label, model_type, [])
    
    disp(['Model ', model_label, ' is completed!! ', progress_count])
end

clear to_disk