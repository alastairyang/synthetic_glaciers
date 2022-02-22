function data_out = re_resol_alongflow(b, h, s, max_len)
%RE_RESOL_ALONGFLOW restore the resolution of straightened geometry along the flow
%   such that the length is 20 km long.
%   Along-flow direction is row, in the data space.
%   This function does not interpolate; it merely swapped out the original
%   coordinate
%
%   Input:
%       s[double, matrix]: surface elevation
%       b[double, matrix]: base elevation
%       max_len[double]:   max cutoff length for flowline in meter
    
    interval = 150; % bedmachine resolution = 150 m
    % let row be x direction, column be y direction
    Nx = size(b, 2);
    Ny = size(b, 1);
    % construct the coordinate
    x = linspace(0, max_len, Nx);
    X = repmat(x, Ny, 1);
    y = transpose(0:interval:(Ny-1)*interval);
    Y = repmat(y, 1, Nx);
    
    % bundle everything into a structure
    data_out.X = X;
    data_out.Y = Y;
    data_out.b = b;
    data_out.h = h;
    data_out.s = s;

end

