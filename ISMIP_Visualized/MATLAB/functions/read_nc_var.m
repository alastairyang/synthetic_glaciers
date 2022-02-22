function all_data = read_nc_var(filepaths, var)
%READ_NC read the .nc at the given path and output requested data
%   filepaths are from search_nc.m
%   outputs include the transformed polar stereographic coordinates (m)
%   as well as the variable of interest var.
%   Now as it is, it only support querying one var
        
    % check how many files to take care of
    N = numel(filepaths); 
    assert(isa(var, 'string') | isa(var, 'char'),...
           'Input var should be a string or char');
    all_data = cell(N, 1);
    
    for i = 1:N
        % save everything for a .nc file as a structure
        filepath = filepaths{i};
        
        t = ncdateread(filepath, 'time');
        Z = ncread(filepath, var);
        inst = ncreadatt(filepath, '/', 'comment');
        source = ncreadatt(filepath, '/', 'source');
        comment = ncreadatt(filepath, '/', 'comment');
        try
            lon = double(ncread(filepath, 'lon'));
            lat = double(ncread(filepath, 'lat'));
            [X, Y] = ll2psn(lat, lon);
            
        % there is a chance that model uses x and y cartesian coords
        % and hence lon, lat do not exist
        catch 
            disp('  ...well there is an error')
            X = double(ncread(filepath, 'x'));
            Y = double(ncread(filepath, 'y'));
            % if it smoothly runs upto here, then it is a coordinate prob.
            
            % we need then to check if X and Y are given as vectors
            if ismember(1, size(X)) && ismember(1, size(Y)) % if so for both
                [X, Y] = meshgrid(X, Y);
                X = X';
                Y = Y';
                disp('   successfully meshgrided')
            end
            disp('   ...Ok polar stereographic proj is given')
        end
        
        % save this structure into one entry of the all_data cell
        ds.X = X;
        ds.Y = Y;
        ds.var = Z;
        ds.lon = lon;
        ds.lat = lat;
        ds.att.institution = inst;
        ds.att.source = source;
        ds.att.comment = comment;
        ds.path = filepath;
        ds.time = t;
        
        all_data{i} = ds;
        
        % print if completed 
        sprintf('file number %d is completed', i)
    end

end

