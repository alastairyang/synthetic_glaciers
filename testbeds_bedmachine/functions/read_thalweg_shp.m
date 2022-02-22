function thalweg_trans = read_thalweg_shp(file_path)
%READ_THALWEG_SHP read in thalweg shapefile and transform the coordinates

    % read in
    thalweg = shaperead(file_path);
    N = size(thalweg, 1);
    
    % Make a nested structure with each fieldname being the glacier id
    for i = 1:N
        [x, y] = ll2psn(thalweg(i).Y, thalweg(i).X);
        thalweg_xy = [x', y'];
        % last rows is NaN, trim out
        thalweg_xy = thalweg_xy(1:end-1,:);
        glaciername = ['glacier','0',num2str(thalweg(i).id)];
        thalweg_trans.(glaciername).('thalweg_xy') = thalweg_xy;
    end
    
    
end

