% This file plot the strainghtened fjord valley
mesh_folder_path = 'Plots/Try/mesh_fig_bh';
straight_pcolor_folder_path = 'Plots/Try/pcolor_png_base';

keyword = 'thalweg';
current_path = pwd;
cd(mesh_folder_path)
all_mesh_thalwegnames = dir(['*',keyword,'*']);

% read in all thalweg files
% get the data from these .fig files
N_thalweg = size(all_mesh_thalwegnames,1);
for i = 1:N_thalweg
    this_name = all_mesh_thalwegnames(i).name;
    openfig(this_name)
    a = get(gca, 'Children');
    xdata = get(a, 'XData');
    ydata = get(a, 'YData');
    zdata = get(a, 'ZData');
    field_name = this_name(1:end-4);
    thalweg_data.(field_name).X = xdata{3};
    thalweg_data.(field_name).Y = ydata{3};
    thalweg_data.(field_name).b = zdata{3};
    thalweg_data.(field_name).s = zdata{2};
end

% use re_resol_alongflow to restore the spatial resolution (at BedMachine)
thalweg_names = fieldnames(thalweg_data);
max_len = 450; % just a big number < 500 to be recognized by the function as a km unit

%% plot the first N
figure
total_plot = 25;
plot_n = 1:total_plot;
plot_dim  = [5,5];

for i = 1:total_plot
    subplot(plot_dim(1), plot_dim(2), i);
    X = thalweg_data.(thalweg_names{plot_n(i)}).X;
    Y = thalweg_data.(thalweg_names{plot_n(i)}).Y;
    b = thalweg_data.(thalweg_names{plot_n(i)}).b;
    %s = thalweg_data_restored.(thalweg_names{i}).s;
    mesh(X,Y,b); 
end