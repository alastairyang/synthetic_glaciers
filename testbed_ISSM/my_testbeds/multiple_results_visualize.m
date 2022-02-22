%% Load the model output data
load('results/sensitivity_3x3.mat')
model_index = 3; % model index number
model_type  = 't_sensitive';
testnames = fieldnames(model_out);

% get the labels for each tests
% load sens data
[geometry, ~, ~, sens_data] = query_data(model_index, model_type);
sens_data = load(sens_data{1});
dataname = ['sens_',num2str(model_index)];
N_test = sens_data.(dataname).N_test;
N_vars  = sens_data.(dataname).N_vars;
varnames = fieldnames(sens_data.(dataname));
titlestrs = [];
for i = 1:N_vars
    testnames = fieldnames(sens_data.(dataname).(varnames{i}));
    for j = 1:N_test
        this_title = varnames{i};
        meanval = mean(sens_data.(dataname).(varnames{i}).(testnames{j}),"all");
        this_title = convertCharsToStrings([this_title, ' = ' , num2str(meanval)]);
        titlestrs = [titlestrs, this_title];
    end
end
    
%% Animate change in lateral profile
[x, y, X, Y,...
 ~, ~, ~, ~,...
 ~, ~,...
 ~, ~,...
 ~, ~, ~, ~, ~, ~, ~] = testbed_data(geometry{1});
if rem(size(X,1), 2) == 0
    mid_i = size(X,1)/2;
else
    mid_i = (size(X,1)+1)/2;
end
thalweg_x = X(mid_i,:);

% construct model testnames
testnames = fieldnames(model_out);
N_test = numel(testnames);
dimension = [2,3];

figure;
for j = 1:N_test
    subplot(dimension(1), dimension(2), j);
    md = model_out.(testnames{j});
    nt = md.timestepping.final_time/md.timestepping.time_step;
    selected = 1:floor(nt*0.05):nt;
    
    % create gradient color palatte
    blueRamp = transpose(linspace(0, 1, numel(selected)));
    black    = zeros(size(blueRamp));
    colors   = [black, black, blueRamp];

    % plot lines for one subplot
    count_i = 0;
    for i = selected
        count_i = count_i + 1;
        thalweg_h_grid = griddata(md.mesh.x, md.mesh.y,...
                                            md.results.TransientSolution(i).Surface,...
                                            X, Y);
        thalweg_base_grid = griddata(md.mesh.x, md.mesh.y,...
                                            md.results.TransientSolution(i).Base,...
                                            X, Y);
        thalweg_h_line = thalweg_h_grid(mid_i,:);
        thalweg_base_line = thalweg_base_grid(mid_i,:);
        thalweg_both = [thalweg_h_line; thalweg_base_line];

        plot(thalweg_x, thalweg_both,'color',colors(count_i,:)); hold on;
    end
    title(titlestrs{j}); 
end

%% Apply EOF to the surface signals
% We collect all data into a 3-D matrix
figure;
for j = 1:N_test
    md = model_out.(testnames{j});
    sampled_t = 1:floor(0.02*nt):nt;
    Ht = zeros(numel(y), numel(x), numel(sampled_t));
    count_i = 0;
    for i = sampled_t
        this_Ht = griddata(md.mesh.x, md.mesh.y,...
                           md.results.TransientSolution(i).Surface,...
                           X, Y);
        count_i = count_i + 1;
        Ht(:,:,count_i) = this_Ht;
    end
    % EOF
    [eofmaps, pc, expvar] = eof(Ht,5);

    % Store the EOF data into nested structures
    EOF.(testnames{j}).eofmaps = eofmaps;
    EOF.(testnames{j}).pc      = pc;
    EOF.(testnames{j}).expvar  = expvar;

    % Reconstruct the timeseries from the first 5 PCs
    t = sampled_t*0.1;
    sum_h = zeros(size(t));
    for i = 1:5
        sum_h = sum_h + pc(i,:)*0.01*expvar(i);
        % store it too
        EOF.(testnames{j}).recons_Ht = sum_h;
    end
    subplot(2,3,j)
    plot(t, sum_h)
    title(titlestrs{j})

end
%% Plot principal components time series and corresponding variabilities
% define colorbar limit
min_map = min(EOF.model_3_sens_test1.eofmaps,[],'all');
max_map = max(EOF.model_3_sens_test1.eofmaps,[],'all');

figure;
subsubplot(5,1,1)
imagesc(x, y, EOF.model_3_sens_test1.eofmaps(:,:,1))
box off
axis tight
ylabel(['pc1: ',num2str(expvar(1)),'%'])
colorbar
caxis([min_map, max_map])

subsubplot(5,1,2)
imagesc(x, y, EOF.model_3_sens_test1.eofmaps(:,:,2))
box off
axis tight
ylabel(['pc2: ',num2str(expvar(2)),'%'])
colorbar
caxis([min_map, max_map])

subsubplot(5,1,3)
imagesc(x, y, EOF.model_3_sens_test1.eofmaps(:,:,3))
box off
axis tight
ylabel(['pc3: ',num2str(expvar(3)),'%'])
colorbar
caxis([min_map, max_map])

subsubplot(5,1,4)
imagesc(x, y, EOF.model_3_sens_test1.eofmaps(:,:,4))
box off
axis tight
ylabel(['pc1: ',num2str(expvar(4)),'%'])
colorbar
caxis([min_map, max_map])

subsubplot(5,1,5)
imagesc(x, y, EOF.model_3_sens_test1.eofmaps(:,:,5))
box off
axis tight
ylabel(['pc1: ',num2str(expvar(5)),'%'])
colorbar
caxis([min_map, max_map])


%% Plot map of principal components
figure;
subsubplot(5,1,1)
plot(t, pc(1,:)/pc(1,1))
box off
axis tight
ylabel(['pc1: ',num2str(expvar(1)),'%'])
title 'The first two PCs'

subsubplot(5,1,2)
plot(t, pc(2,:)/pc(2,1))
box off
axis tight
set(gca,'yaxislocation','right')
ylabel(['pc2: ',num2str(expvar(2)),'%'])

subsubplot(5,1,3)
plot(t, pc(3,:)/pc(3,1))
box off
axis tight
ylabel(['pc3: ',num2str(expvar(3)),'%'])

subsubplot(5,1,4)
plot(t, pc(4,:)/pc(4,1))
box off
axis tight
set(gca,'yaxislocation','right')
ylabel(['pc4: ',num2str(expvar(4)),'%'])

subsubplot(5,1,5)
plot(t, pc(5,:)/pc(5,1))
box off
axis tight
ylabel(['pc5: ',num2str(expvar(5)),'%'])

%% plot time series at all points
coor_N = numel(x)*numel(y);
t = 0:0.1:(nt-1)*0.1;
all_timeseries = zeros(coor_N, nt);
figure
for i = 1:coor_N
    % get the corresponding subscript index
    [row, col] = ind2sub([numel(y), numel(x)], i);
    this_timeserie = Ht(row, col,:);
    this_timeserie = this_timeserie(:);
    this_timeserie = this_timeserie/this_timeserie(1);
    all_timeseries(i,:) = this_timeserie;
    plot(t, this_timeserie)
    hold on;
end

    