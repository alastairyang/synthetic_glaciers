function map_masked = polygon2mask(map, vertices)
%POLYGON2SQUAREMASK create a rectangular mask 
%
%   Input:
%       map: a structure including X, Y, and other glacier variables "var_i"
%       vertices: a Nx2 matrix where first column is x coordinate and the
%           second column is y coordinate
%
%   Output:
%       map_mask: a structure including X, Y, and other variables of the
%           region of interest.

    % First confirm that there are X and Y
    map_names = fieldnames(map);
    assert(sum(cellfun(@(x) x=='X', map_names)), 'There is no X!');
    assert(sum(cellfun(@(x) x=='Y', map_names)), 'There is no Y!');
    % Then assign X, Y
    X = map.X;
    Y = map.Y;
    xi = vertices(:,1);
    yi = vertices(:,2);
    
    size_x = size(X, 2);
    size_y = size(X, 1);
    mask = poly2mask(xi, yi, size_y, size_x);
    
    % sum vertical and horizontal dimension
    % the first and last non-zero entries are the bounds
    % we record their subscript indices
    
    % sum along the first dimension: get the bounds for length along the
    % horizontal direction
    mask_sum_1 = sum(mask, 1);
    % sum along the second dimension: ...along the vertical direction
    mask_sum_2 = sum(mask, 2);
    
    % find the indices of four vertices
    left_idx  = find(mask_sum_1>0, 1, 'first');
    right_idx = find(mask_sum_1>0, 1, 'last');
    up_idx    = find(mask_sum_2>0, 1, 'first');
    down_idx  = find(mask_sum_2>0, 1, 'last');
    
    % we need to crop out both the map (coordinates and vars)
    % and the mask for later use
    
    % Then We will perform element wise multiplication to 0 out the entries not
    % in the original polygon, but we first need to manually make 0s in the
    % map to be non-zero (e.g., 0.00001)
    % then make all zeros NaNs
    sub_0 = 0.0001;
    mask_crop = mask(up_idx:down_idx, left_idx:right_idx);
    
    X_crop = X(up_idx:down_idx, left_idx:right_idx);
    Y_crop = Y(up_idx:down_idx, left_idx:right_idx);
    
    map_masked.X = X_crop;
    map_masked.Y = Y_crop;
    
    % iterate over other field variables
    for i = 1:numel(map_names)
        % skip X and Y as we do not perform masking on them
        if map_names{i} == 'X' || map_names{i} == 'Y'
            continue;
        end
        % let z be a generic variable name
        z = map.(map_names{i});
        z_crop = z(up_idx:down_idx, left_idx:right_idx);
        z_crop(z_crop==0) = sub_0;
        z_masked = mask_crop.*z_crop;
        z_masked(z_masked==0) = NaN;
        % assign to the output structure
        map_masked.(map_names{i}) = z_masked;
    end
    
    
    
end

