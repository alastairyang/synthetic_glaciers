%% Analyze Michael Wood et al (2021) data (no truncation)
path = "/Users/donglaiyang/Library/Mobile Documents/com~apple~CloudDocs/Documents/Documents_Alastair’s_MacBook_Pro/Buffalo/Research/Data and SI/Wood et al., 2021/Original_Data_removeNa.xlsx";
W_data_m = table();
for i = 2:8
    W_data = rows2vars(readtable(path, 'Sheet',i,'ReadRowNames',true));
    W_data_m = [W_data_m; W_data];
end
units = W_data_m(1,:);
W_data_m(1,:) = [];

% convert cells to double matrix
values = W_data_m(:,2:end);
mean_depth = mean(values.MeanDepth, 'all');
mean_fjord_width = mean(values.MeanFjordWidth, 'all');
min_depth = min(values.MeanDepth, [],'all');
min_fjord_width = min(values.MeanFjordWidth, [], 'all');
max_depth = max(values.MeanDepth, [], 'all');
max_fjord_width = max(values.MeanFjordWidth, [],'all');
width_qtl = quantile(values.MeanFjordWidth, [0.025, 0.25 0.50 0.75 0.975])
depth_qtl = quantile(values.MeanDepth, [0.025, 0.25 0.50 0.75 0.975])

figure;
subplot(1,2,1)
histogram(values.MeanDepth)
title('Mean Grounding line depth distribution')
xlabel('m')
subplot(1,2,2)
histogram(values.MeanFjordWidth(values.MeanFjordWidth<20))
title('Mean Fjord width distribution')
xlabel('km')
