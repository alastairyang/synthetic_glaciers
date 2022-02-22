function [] = plot_bed()
%PLOT_H Plot ice thickness from BedMachine V4

    filename = 'bedtopo_GrIS.mat';
    load(filename)
    
    X = bedtopo.X;
    Y = bedtopo.Y;
    B = bedtopo.B;
    
    X_max = max(X(:));
    X_min = min(X(:));
    Y_max = max(Y(:));
    Y_min = min(Y(:));

    ratio = (Y_max - Y_min)/(X_max - X_min);
    fig_width = 1200;
    
    figure('position',[100,100,fig_width, fig_width*ratio]);
    greenland()
    hold on
    mesh(X, Y, B);
    hold off
end

