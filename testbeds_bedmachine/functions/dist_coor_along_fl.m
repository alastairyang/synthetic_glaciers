function [coor_at_dist, index] = dist_coor_along_fl(fl_data, dist)
%DIST_COOR_ALONG_FL for a given list of flow line coordinates and a maximum
%distance, we will output the coordiante that the given distance is reached
%
%   Input:
%       fl_data[array]: flow line coordinate x and y
%       dist[double]:   a single number in km
%
%   Output:
%       coor_at_dist[double]: coordinate at the given distance

    i = 1;
    dist_m = dist*1000;
    dist_computed = 0;
    while dist_computed < dist_m
        if i == 1
            last_data = fl_data(i,:);
            i = i+1;
            continue
        end
        temp_dist = sqrt((fl_data(i,1)-last_data(1,1))^2 + ...
                         (fl_data(i,2)-last_data(1,2))^2);
        dist_computed = dist_computed + temp_dist;
        last_data = fl_data(i,:);
        i = i+1;
        
        % get to the end of array before reaching the required distance
        % just break; deal with whatever we have
        if i+1 > size(fl_data,1)
            break;
        end
    end
    
    % once we are out of this while loop
    % i - 1 is the coordinate we are lookking for
    index = i-1;
    coor_at_dist = fl_data(index,:);
    
end

