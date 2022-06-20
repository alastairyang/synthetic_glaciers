%% This script get data from a synthetic glacier and visualize the forcing, 
% both their spatial extent and the changes in time.
index = '0';
model_type = 't';
[geometry, ~] = query_data(index, model_type);
syn = testbed_data(geometry{1});
x = syn.X(1,:);
y = syn.Y(:,1);
min_fric = min(syn.transient_fric_coef.data{2},[],'all');
max_fric = max(syn.transient_fric_coef.data{1},[],'all');
min_rB   = min(syn.transient_rheoB.data{2},[],'all');
max_rB   = min(syn.transient_rheoB.data{1},[],'all');
%%
figure('Position',[100,100,200,300]);
this_fric = syn.transient_fric_coef.data{2};
imagesc(y/1e3, x/1e3, this_fric')
colorbar
caxis([min_fric, max_fric])
print(gcf,'Graphs/fric_max.png','-dpng','-r300');  

figure('Position',[100,100,200,300]);
this_rB = syn.transient_rheoB.data{2};
imagesc(y/1e3, x/1e3, this_rB')
colorbar
caxis([min_rB, max_rB])
print(gcf,'Graphs/rB_max.png','-dpng','-r300');  

figure('Position',[100,100,300,200]);
melt1 = mean(syn.shelf_melt.transient_melt{1},'all');
fric2 = mean(syn.shelf_melt.transient_melt{2},'all');
melt3 = mean(syn.shelf_melt.transient_melt{3},'all');
t_cell = syn.shelf_melt.melt_years;
taxis = [0, 10, 20, 30, 40, 50, t_cell{1}, t_cell{2}, t_cell{3}, 100];
melts = [melt1,melt1, melt1, melt1, melt1, melt1, melt1, fric2, melt3, melt3];
new_t = 0:10:100;
fric_interp = interp1(taxis, melts, new_t);
plot(new_t, fric_interp,'-*')
xlabel('Time (year)','FontSize',13, 'FontName','Times')
ylabel('melt rate (m/y)','FontSize',13, 'FontName','Times')

print(gcf,'Graphs/meltrate.png','-dpng','-r300');  

%% Time series of magnitude, fric. coef and rheoB
figure('Position',[100,100,300,200]);
fric1 = mean(syn.transient_fric_coef.data{1},'all');
fric2 = mean(syn.transient_fric_coef.data{2},'all');
fric3 = mean(syn.transient_fric_coef.data{3},'all');
t_cell = syn.transient_fric_coef.years;
taxis = [0, 10, 20, 30, 40, 50, t_cell{1}, t_cell{2}, t_cell{3}, 100];
frics = [fric1,fric1, fric1, fric1, fric1, fric1, fric1, fric2, fric3, fric1];
new_t = 0:10:100;
fric_interp = interp1(taxis, frics, new_t);
plot(new_t, fric_interp,'-*')
xlabel('Time (year)','FontSize',13, 'FontName','Times')
ylabel('max Fric.Coef','FontSize',13, 'FontName','Times')
print(gcf, 'Graphs/fric_t.png','-dpng','-r300')

figure('Position',[100,100,300,200]);
fric1 = mean(syn.transient_rheoB.data{1},'all');
fric2 = mean(syn.transient_rheoB.data{2},'all');
fric3 = mean(syn.transient_rheoB.data{3},'all');
t_cell = syn.transient_rheoB.years;
taxis = [0, 10, 20, 30, 40, 50, t_cell{1}, t_cell{2}, t_cell{3}, 100];
frics = [fric1,fric1, fric1, fric1, fric1, fric1, fric1, fric2, fric3, fric1];
new_t = 0:10:100;
fric_interp = interp1(taxis, frics, new_t);
plot(new_t, fric_interp,'-*')
xlabel('Time (year)','FontSize',13, 'FontName','Times')
ylabel('max Rheology B','FontSize',13, 'FontName','Times')
print(gcf, 'Graphs/rheoB_t.png','-dpng','-r300')

