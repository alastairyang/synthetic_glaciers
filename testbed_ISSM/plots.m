% visualizing results
time = cell2mat({md.results.TransientSolution(:).time});
thickness = cell2mat({md.results.TransientSolution(:).Thickness});
thickness_mean = mean(thickness, 1);
vel = cell2mat({md.results.TransientSolution(:).Vel});
vel_mean = mean(vel, 1);

figure;
subplot(1,2,1)
plot(time, thickness_mean, 'DisplayName', 'Thickness')
legend

subplot(1,2,2)
plot(time, vel_mean, 'DisplayName', 'Velocity')
legend