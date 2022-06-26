function [Xq, Yq, Zq] = meshgrid_downsample(X, Y, Z, model_type)
%MESHGRID2OUTLINE downsampling the given meshgrid and variables. This
%function also first informs users the original meshgrid dimensions and
%then ask for the desired new downsampled dimension
%
%   example: [Xq, Yq, Zq] = meshgrid2outline(X, Y, thickness);
%
%   input:
%       X: a meshgrid of X
%       Y: a meshgrid of Y
%       Z: a meshgrid of variable Z
%
%   output:
%       Xq: a downsampled meshgrid of X
%       Yq: a downsampled meshgrid of Y
%       Zq: a downsampled meshgrid of variable Z

    % test dimension agreement
    assert((size(X,1) == size(Y,1)) && (size(X,2) == size(Y,2)),...
            'The dimensions of X,Y do not agree')
    % test X is invariant along column, Y is invariant along row
    if range(X(:,1)) ~= 0
        X = X';
        assert(range(X(:,1))==0, 'X is not invariant along any dimension')
    end
    if range(Y(1,:)) ~= 0
        Y = Y';
        assert(range(Y(1,:))==0, 'Y is not invariant along any dimension')
    end
        
    % requires user's input
%     sprintf('The dimension of the meshgrid: x = %d, y = %d', size(X,2), size(X,1))
%     prompt = 'Downsample to what percent? In decimal. \n';
    if strcmp(model_type, 'ss')
        p = 0.2;
        disp('    The meshgrid is downscaled to 20%')
    elseif strcmp(model_type, 'spinup')
        p = 0.1;
        disp('    The meshgrid is downscaled to 10%')
    end
    nx = floor(size(X,2)*p);
    ny = floor(size(X,1)*p);
    
    % get the min and max values for X and Y
    max_X = max(X(1,:));
    min_X = min(X(1,:));
    max_Y = max(Y(:,1));
    min_Y = min(Y(:,1));
    % test again that min and max from one vector should also be global
    % max and min
    assert(max_X == max(X(:)), 'X max values do not match.')
    assert(max_Y == max(Y(:)), 'Y max values do not match.')
    
    Xq = linspace(min_X, max_X, nx);
    Yq = linspace(min_Y, max_Y, ny);
    [Xq, Yq] = meshgrid(Xq, Yq);
    Zq = interp2(X, Y, Z, Xq, Yq);
    
end

