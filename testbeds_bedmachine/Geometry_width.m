
%% Get width distribution and plan view topography data for visualization
widths_26 = width_distribution('glacier026_thalweg_xy_mesh_bh.fig');
widths_1  = width_distribution('glacier01_thalweg_xy_mesh_bh.fig');
widths_8  = width_distribution('glacier08_thalweg_xy_mesh_bh.fig');

% plot the plan view and the distribution of width
bed_26 = get_bed('glacier026_thalweg_xy_mesh_bh.fig');
bed_1  = get_bed('glacier01_thalweg_xy_mesh_bh.fig');
bed_8  = get_bed('glacier08_thalweg_xy_mesh_bh.fig');

%% Import images
image_26 = imread('Plots/Try/pcolor_png_original_base/glacier026_thalweg_xy_pcolor_original_b.png');
image_1  = imread('Plots/Try/pcolor_png_original_base/glacier01_thalweg_xy_pcolor_original_b.png');
%% make grid plots
figure;
% % find width limits
% subplot(2,2,1)
% imshow(image_1)
% title('Glacier # 1','FontSize',12)
% xlabel('x (m)')
% ylabel('y (m)')
% 
% subplot(2,2,2)
% imshow(image_26)
% title('Glacier # 26','FontSize',12)
% xlabel('x (m)')
% ylabel('y (m)')

%subplot(2,2,[3,4])
histogram(widths_1(widths_1>0)); hold on;
histogram(widths_26(widths_26>0)); hold off;
xlim([min([min(widths_26), min(widths_1)]), max([max(widths_26), max(widths_1)])])
xlabel('Width (meter)','FontSize',13, 'FontName','Times')
ylabel('Frequency','FontSize',13, 'FontName','Times')
leg1 = legend('Sermeq Kujalleq','Upernavik Isstrøm','FontSize',13, 'FontName','Times');
set(leg1,'Box','off')

print(gcf,'Graphs/two_glacier_widths.png','-dpng','-r300');  

% subplot(2,4,5)
% histogram(widths_26(widths_26>0))
% xlim([min([min(widths_26), min(widths_1)]), max([max(widths_8), max(widths_1)])])
% xlabel('Width (meter)')
% ylabel('Frequency')

% subplot(2,3,3)
% h = pcolor(bed_8.X, bed_8.Y, bed_8.b);
% set(h, 'EdgeColor','none')
% title('Glacier # 8','FontSize',12)
% xlabel('x (m)')
% ylabel('y (m)')
% 
% subplot(2,3,6)
% histogram(widths_8(widths_8>0))
% xlim([min([min(widths_8), min(widths_1)]), max([max(widths_8), max(widths_1)])])
% xlabel('Width (meter)')
% ylabel('Frequency')

%% APPENDIX: Functions
function widths = width_distribution(filename)
%%WIDTH_DISTRIBUTION Estimate the width distribution from straightened
%%glaciers (across-flow cross-sectional profiles)
    fig_path = ['Plots/Try/mesh_fig_bh/', filename];
    openfig(fig_path);

    % get data from figure
    a = get(gca, 'Children');
    xdata = get(a, 'XData');
    ydata = get(a, 'YData');
    zdata = get(a, 'ZData');

    X = xdata{3};
    Y = ydata{3};
    b = zdata{3};

    % Use Findpeaks from signal processing toolbox
    dy = 150; 
    y = Y(:,1);
    N_col = size(b, 2);
    N_row = size(b, 1);
    max_w = zeros(N_col,1);

    y_center = y((N_row+1)/2,1);

    for i = 1:100
        % uncomment the following line to visualize the profile
        % findpeaks(-1*b(:,i),y,'Annotate','extents','WidthReference','halfprom');
        [pks,locs,w,~] = findpeaks(-1*b(:,i),y,'Annotate','extents','WidthReference','halfprom');
        if isempty(pks)
            continue
        end
        % Three scoring criterions: 
        % 1. largest peak height (deepest trough)
        % 2. Closest to the center line
        % 3. Largest width
        [~, I_pks]   = sort(pks, 'descend');
        [~, I_locs] = sort((locs-y_center).^2, 'ascend');
        [~, I_w]       = sort(w, 'descend');
        % we only consider the first 3
        if length(I_pks)<3 || length(I_locs)<3 || length(I_w)<3
            % choose the shortest length
            n = min([length(I_pks), length(I_locs), length(I_w)]);
        else
            n = 3;
        end
        all_I   = [I_pks(1:n)', I_locs(1:n)', I_w(1:n)'];
        all_I_u = unique(all_I);
        credits = repmat(n:-1:1,1,3);
        scores  = zeros(size(all_I_u));
        for j = 1:length(all_I)
            idx = find(all_I_u == all_I(j));
            scores(idx) = scores(idx) + credits(j);
        end
        [~, I_max] = max(scores);
        max_w(i) = w((all_I_u(I_max)));
    end
    
    widths = max_w;
end

function bed = get_bed(filename)
%%GET_BED get the bed data from the figure file
    fig_path = ['Plots/Try/mesh_fig_bh/', filename];
    openfig(fig_path);

    % get data from figure
    a = get(gca, 'Children');
    xdata = get(a, 'XData');
    ydata = get(a, 'YData');
    zdata = get(a, 'ZData');

    X = xdata{3};
    Y = ydata{3};
    b = zdata{3};
    
    bed.X = X;
    bed.Y = Y;
    bed.b = b;

end

% %% FFT method
% dy = 150; 
% y = Y(:,1);
% N_col = size(b, 2);
% 
% signal_len = size(b,1);
% PSD = zeros(size(b));
% n = size(b,1);
% fft_out = fft(b,n);
% for i = 1:N_col
%     %%PSD(:,i) = fft_out(:,i).*conj(fft_out(:,i))/n;
%     PSD(:,i) = abs(fft_out(:,i));
% end
% figure;
% freq = (max(y)-min(y))/(dy*n)*(1:n);
% perd = 1*(max(y)-min(y))./freq/2;
% plot(perd, PSD)
% 
% % find the second largest peak (second largest number)
% % period also > 100
% PSD_trim = PSD(perd>100,:);
% [B, I] = maxk(PSD_trim, 2, 1);
% B_2nd = B(2,:); 
% I_2nd = I(2,:);
% perd_2nd = perd(I_2nd);
% B_1st = B(1,:);
% I_1st = I(1,:);
% perd_1st = perd(I_1st);
% 
% figure; histogram(perd_2nd); hold on; histogram(perd_1st)
% legend('2nd max','1st max')
% xlabel('Width (m)')
% ylabel('Count')
% title('Max period count')
