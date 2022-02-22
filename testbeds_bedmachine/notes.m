%% open and visualize fig
openfig('glacier0003_flowline03_mesh_bh.fig')
set(gcf, 'visible','on')

%% get data from figure
a = get(gca, 'Children');
xdata = get(a, 'XData');
ydata = get(a, 'YData');
zdata = get(a, 'ZData');

%% plot
figure; mesh(xdata{1}, ydata{1}, zdata{1});
figure; mesh(xdata{2}, ydata{2}, zdata{2});
figure; mesh(xdata{3}, ydata{3}, zdata{3});

