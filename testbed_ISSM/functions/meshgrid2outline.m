function [] = meshgrid2outline(X, Y)
%%This function outputs outline coordinates from a meshgrid. The input
%%coordinate system should be cartesion, meaning the bottom-left element
%%should be (0,0)
%
%   example: [outline_coor, total_N] = meshgrid2outline(X, Y);
%
%   input:
%       X: a meshgrid of X
%       Y: a meshgrid of Y
%
%   output:
%       No functional output, but saves a text file.
%   
%   We define the total number of sampled points to be (nx-1)*2 + (ny-1)*2

    N_x = size(X, 2); % number of elements of each row
    N_y = size(Y, 1); % number of elements of each column
    % first we check if the coordinate system is cartesian
    % we only need to check y mesh since x is the same regardless of
    % coordinate system
    if Y(N_y, 1) ~= 0
        % then we flip Y mesh
        Y = flip(Y, 1);
    end
    
    %%%%%%%
    % We keep the original vertices coordinates
    % bottom ind_X = 1 to N_x-1
    bottom_x = X(1, 1:N_x-1); % X is homogenous along the first dim
    bottom_y = Y(N_y, 1:N_x-1);
    bottom_xy = [bottom_x', bottom_y'];
    
    
    % right side, index_X = N_x, index_Y = 2 to N_y (skip the top right one)
    right_y = Y(2:N_y, N_x);
    right_x = X(1, N_x)*ones(numel(right_y), 1);
    % reverse order so that it is from bottom up
    right_xy = [right_x, right_y];
    right_xy_flip = flip(right_xy, 1);
    
    
    % top side, similar to the bottom, except starting at the second to the
    % second last
    top_x = X(1, 2:N_x);
    top_y = Y(1, 2:N_x);
    % reverse order so that it is from right to left
    top_xy = [top_x', top_y'];
    top_xy_flip = flip(top_xy, 1);
    
    
    % left side...
    left_y = Y(1:N_y, 1);
    left_x = X(1, 1)*ones(numel(left_y), 1);
    left_xy = [left_x, left_y];
    
    % append all arrays
    outline_coor = [bottom_xy; right_xy_flip; top_xy_flip; left_xy];
    
    % test if the last entry row is (0,0) which is required
    % if false, add a row of zeros
    if ~isequal(outline_coor(1,:), outline_coor(end, :))
        outline_coor = cat(1, outline_coor, outline_coor(1,:));
    end
    
    % calculate total number of entries
    total_N = size(outline_coor, 1);
    
    %% write a txt file
    fileID = fopen('Domain.exp','w');
    fprintf(fileID,'%s\n','## Name:DomainOutline');
    fprintf(fileID,'%s\n','## Icon:0');
    fprintf(fileID,'%s\n','# Points Count  Value');
    fprintf(fileID,'%f %f\n', total_N, 1);
    fprintf(fileID,'%s\n','# X pos Y pos');
    fprintf(fileID,'%f %f\n', outline_coor.');
    fclose(fileID);
    
    %writematrix(outline_coor, 'Domain.txt', 'delimiter', ' ','precision','%.4f');
end

