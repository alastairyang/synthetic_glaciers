function [centroids_x, centroids_y] = triangle_centroid(element_number, x, y)
%TRIANGLE_CENTROID This function calculates the centroid location of the
%triangular elements
%
%   Input:
%       element_number [double array]: Mx3 matrix recording the index number
%                                     of x,y coordianates of each vertex
%       x              [double array]: Nx1 vector recording the x coor
%       y              [double array]: Nx1 vector recording the y coor
%       NOTE: In general, M > N
%   Output:
%       centroids      [double array]: Mx2 matrix recording the centrod
%                                      locations (x,y coor).

    N_element = size(element_number, 1);
    centroids = zeros(size(N_element,1),2);
    for i = 1:N_element
        vertex1 = [x(element_number(i,1)), y(element_number(i,1))];
        vertex2 = [x(element_number(i,2)), y(element_number(i,2))];
        vertex3 = [x(element_number(i,3)), y(element_number(i,3))];
        vertices = [vertex1; vertex2; vertex3];
        % from means of x, y -> centroid
        centroids(i,:) = mean(vertices, 1);
    end
    centroids_x = centroids(:,1);
    centroids_y = centroids(:,2);
end

