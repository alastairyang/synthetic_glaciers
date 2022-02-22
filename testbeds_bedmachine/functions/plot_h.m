function [] = plot_h()
%PLOT_H Plot ice thickness from BedMachine V4

    filename = 'thickness_GrIS.mat';
    load(filename)
    
    X = thickness.X;
    Y = thickness.Y;
    h = thickness.h;
    
    X_max = max(X(:));
    X_min = min(X(:));
    Y_max = max(Y(:));
    Y_min = min(Y(:));

    ratio = (Y_max - Y_min)/(X_max - X_min);
    fig_width = 1200;
    
    figure('position',[100,100,fig_width, fig_width*ratio]);
    greenland()
    hold on
    mesh(X, Y, h);
    hold off
end

