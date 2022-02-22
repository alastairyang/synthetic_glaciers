function bm_data = import_bm(dir)
%import_bd Import bedmachine data from a given directory
%   
%   example:  import_bd('my_directory')
%
%   Input:
%       dir: directory to the file. A string or characters
%
%   Output:
%       bm_data: structure containing relevant bedmachine data

    x = double(ncread(dir, 'x'));
    y = double(ncread(dir, 'y'));
    [X, Y] = meshgrid(x, y);
    s = double(ncread(dir, 'surface'));
    h = double(ncread(dir, 'thickness'));
    b = double(ncread(dir, 'bed'));
    mk = double(ncread(dir, 'mask'));
    
    % into the structure
    bm_data.X = X;
    bm_data.x = x;
    bm_data.Y = Y;
    bm_data.y = y;
    bm_data.s = s';
    bm_data.h = h';
    bm_data.b = b';
    bm_data.mk = mk';
    
end

