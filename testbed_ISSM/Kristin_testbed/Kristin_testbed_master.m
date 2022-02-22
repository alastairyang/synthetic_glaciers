% This file runs multiple ISMIP tests that differ from each other slightly
clear
clc
clear params

tic
md = model_execute();
view(45, 10)
    
% plot and save figure
plotmodel(md,'data', md.results.StressbalanceSolution.Vel, 'colorbar','on','edgecolor','w');
view(45, 10)
