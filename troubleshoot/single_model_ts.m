% Single Model Troubleshoot (ts)
% Problem: some models have wiggles at the glacier bed
% here we are using syn_00, since it has the most wiggles

geom_path = '/Users/donglaiyang/Library/Mobile Documents/com~apple~CloudDocs/Documents/Documents_Alastair’s_MacBook_Pro/Buffalo/Research/git_research/testbeds_bedmachine/Synthetic glaciers/syn_00.mat';
model_index = '00';

%% Steady state
ts_my_model_execute_ss(geom_path, model_index, 'ss');

%% Spin up
% velocity data from steady state is stored locally
model_type = 'spinup';
vel_path = ['steady state/ss_V_',model_index,'.mat'];
md_spinup = ts_my_model_execute_t(geom_path, vel_path, [], model_index, model_type, []);