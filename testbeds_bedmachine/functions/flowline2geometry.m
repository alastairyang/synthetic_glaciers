function data_all = flowline2geometry(fl_data, map_data, data_res, sub_fl, method)
%FLOWLINE2BED extract geometry (bed, ice thickness) from BedMachine using
%the specified method
%
%   Input:
%       fl_data[struc]:  flowline dataset.
%       map_data[struc]: BedMachine dataset.
%           map_data.X
%           map_data.Y
%           map_data.b
%           map_data.h
%           map_data.s
%           map_data.mk
%       data_res[char]: 
%           'high'    : flow line data has higher resolution; needs to be
%                       coarsen
%           'low'     : flow line (or thalweg here) data has lower resolution;
%                       needs to be refined
%       sub_fl[double array]: 
%       method[char]:    method for extraction
%           'original': keep the curved flowline and extract the minimal
%               boundary rectangle
%           'straight': use coordinate transformation and interpolation to
%               produce a straight glacier channel
%
%   Output:
%       data_all[struc]: nested struc storing geometry info for each flowline

    %% We sample one flow line per glacier
    if ~strcmp(method, 'straight')
        error('Currently we do not support this method')
    end
    
    glacier_names = fieldnames(fl_data);
    N_glacier = numel(glacier_names);
    for i = 1:N_glacier
        % number of flowlines
        flowline_names = fieldnames(fl_data.(glacier_names{i}));
        N_fl = numel(flowline_names);
        % choose the subflowline
        %midnum_fl = middle_number(N_fl);
        
        % if using my own handpicked thalweg, there is always only one line
        if strcmp(data_res, 'low')
            subfl = 1;
        else
            subfl = sub_fl;
        end
        % combine the glacier name and flowline name
        % and save to a new structure
        sampled_fl_data.([glacier_names{i}, '_', flowline_names{subfl}])...
                  = fl_data.(glacier_names{i}).(flowline_names{subfl});
    end
    
    % we obtain a structure of sampled flowline data
    
    %% Extract the geometry using the sampled flow line data
    ds = 2000; % max flow line length in km for our extraction
    n_rect = 150;
    n_sqr  = n_rect/2; 
    
    % unpack map_data
    X = map_data.X;
    Y = map_data.Y;
    b = map_data.b;
    h = map_data.h;
    s = map_data.s;
    mk = map_data.mk;
    
    % EXTRACTION!
    % start iteration for all sampled flowlines
    sp_flowline_names = fieldnames(sampled_fl_data);
    N_sp_flnames = numel(sp_flowline_names);
    for i = 1:N_sp_flnames
        tic
        
        this_fl = sampled_fl_data.(sp_flowline_names{i});
        % first, crop out a small part of BedMachine, otherwise a gird search
        % of points will take way too long
        % cc stands for crop_crop, which covers the flowline limited by ds
        [X_c, Y_c, b_c, h_c, s_c, mk_c,...
         X_cc, Y_cc, b_cc, h_cc, s_cc, mk_cc] =...
                       rect_from_point2(this_fl, X,Y,b,h,s, mk, n_rect, ds);
        
        if strcmp(data_res, 'high') % Dennis's flow line data
            % then we need to coarsen the flowline resolution since it is
            % higher than the BedMachine
            fl_data_coarse = coarsen_fl_points(this_fl, X_c, Y_c);
            fl_data_ongrid_x = fl_data_coarse(:,1);
            fl_data_ongrid_y = fl_data_coarse(:,2);
        else % my thalweg data
            fl_data_refine = refine_fl_points(this_fl, X_c, Y_c);
            fl_data_ongrid_x = fl_data_refine(:,1);
            fl_data_ongrid_y = fl_data_refine(:,2);
        end
        
        % finally, we invoke straighten_flowline function
        % to transform the curvilinear coordinate of flowline
        [X_out, Y_out, b_out, h_out, s_out] = straighten_flowline(fl_data_ongrid_x, fl_data_ongrid_y,...
                                X_c, Y_c, b_c, h_c, s_c, mk_c,...
                                n_sqr, ds);
                            
        %% make and save mesh plots
        m = figure('visible','off');
        mesh(X_out,Y_out,b_out); hold on; % base
        mesh(X_out, Y_out, s_out); hold on; % surface
        mesh(X_out, Y_out, s_out - h_out - b_out); hold off; % zero in land ice; non-zero in marine-terminating section
        saveas(m, ['Plots/Try/mesh_fig_bh/',sp_flowline_names{i},'_mesh_bh.fig'])
        
        % make and save pcolor plots after transform., only for the base though
        p = figure('visible','off');
        a = pcolor(X_out, Y_out, b_out);
        set(a, 'EdgeColor', 'none'); hold off
        saveas(p, ['Plots/Try/pcolor_png_base/',sp_flowline_names{i},'_pcolor_b.png'])
        
        % make and save pcolor plots of the original, 
        % with dot at max length marked with a different color
        coor_max_dist = dist_coor_along_fl(this_fl, ds);
        q = figure('visible','off');
        qp = pcolor(X_cc, Y_cc, b_cc); set(qp, 'EdgeColor','none'); hold on
        scatter(coor_max_dist(1,1), coor_max_dist(1,2), 15, 'b','filled'); hold on
        scatter(fl_data_ongrid_x, fl_data_ongrid_y, 3, 'r','filled'); hold off
        saveas(q, ['Plots/Try/pcolor_png_original_base/',sp_flowline_names{i},...
                   '_pcolor_original_b.png']);
        
        %% save output into structures
        % cropped, non-transformed
        data_crop.X = X_cc;
        data_crop.Y = Y_cc;
        data_crop.b = b_cc;
        data_crop.h = h_cc;
        data_crop.s = s_cc;
        data_crop.mk = mk_cc;
        
        % cropped, and transformed
        data_trans.X = X_out;
        data_trans.Y = Y_out;
        data_trans.b = b_out;
        data_trans.h = h_out;
        data_trans.s = s_out;
        
        data_all.(sp_flowline_names{i}).('crop') = data_crop;
        data_all.(sp_flowline_names{i}).('trans') = data_trans;
        
        % print that this flow line is completed
        disp([sp_flowline_names{i}, ' is now complete!'])
        
        toc
    end
end

