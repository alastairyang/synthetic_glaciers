function syn = testbed_data(directory)
%TESTBED_DATA This function rename the synthetic glacier data structure
%   Input:
%       directory[string/char]: absolute path to the geometry file
%
%   Ouput:
%       syn[structure]: same data (add x,y) but renames (does not carry the
%                       specific glacier index.
    
    data = load(directory);
    name = fieldnames(data);
    syn = data.(name{1});
    
%     This function only accepts meshgrid coordinates
%     try
%         size_x = size(syn.x);
%         if ismember(1,size_x) % it is a vector
%             disp('Only meshgrid coordinates are accepeted')
%             quit
%         end
%     catch
%         size_x = size(syn.X);
%         if ismember(1,size_x) % it is a vector
%             disp('Only meshgrid coordinates are accepeted')
%             quit
%         end
%     end
%     
%     parse the data and convert from km to m
%     print mean of each value to make sure that the magnitudes are
%     consistent for all variables
%     if max(syn.s,[],'all')<100 ||...
%        max(syn.bed,[],'all')<100 ||...
%        max(syn.h,[],'all')<100
%         s = double(syn.s*1000.0);
%         sprintf('Mean surface elevation is %f meter', mean(s, 'all'))
%         bed = double(syn.bed*1000.0);
%         sprintf('Mean basal elevation is %f meter', mean(bed, 'all'))
%         H = double(syn.h*1000.0);
%         sprintf('Mean ice thickness is %f meter', mean(H, 'all'))
%     else
%         s = double(syn.s);
%         sprintf('Mean surface elevation is %f meter', mean(s, 'all'))
%         bed = double(syn.bed);
%         sprintf('Mean basal elevation is %f meter', mean(bed, 'all'))
%         H = double(syn.h);
%         sprintf('Mean ice thickness is %f meter', mean(H, 'all'))
%     end
%     
%     detect if units are km or m
%     assume that value < 100, it is km
%     if max(syn.Y,[],'all') < 100 || max(syn.X,[],'all') < 100
%         Y_g = double(syn.Y*1000.0);
%         sprintf('Y min is %f, and max is %f', min(Y_g(:,1)), max(Y_g(:,1)))
%         X_g = double(syn.X*1000.0);
%         sprintf('X min is %f, and max is %f', min(X_g(1,:)), max(X_g(1,:)))
%     else
%         Y_g = double(syn.Y);
%         sprintf('Y min is %f, and max is %f', min(Y_g(:,1)), max(Y_g(:,1)))
%         X_g = double(syn.X);
%         sprintf('X min is %f, and max is %f', min(X_g(1,:)), max(X_g(1,:)))
%     end

        
    % get the x,y coordinate vector
    % x is along flow-direction and is row
    x = syn.X(1,:);
    y = syn.Y(:,1);
    
    syn.x = x;
    syn.y = y;


end