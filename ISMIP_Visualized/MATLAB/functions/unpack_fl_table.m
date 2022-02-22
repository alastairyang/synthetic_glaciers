function [coor, names, unique_names] = unpack_fl_table(fl_table)
%UNPACK_FL_TABLE Unpacking flow line data table into double and string
%arrays
%
%   example: [coor, names, unique_names] = unpack_fl_table(fl_table)
%  
%   Input:
%       fl_table: flow line data table, consisting of X, Y, and names
%
%   Output:
%       coor: Nx2 array cartesian coordinates
%       names: glacier label string arrays
%       unique_names: the unique glacier labels in the names string array

    assert(size(fl_table,2) == 3, 'The number of table column is not three!')
    
    unique_names = unique(fl_table{:, 3});
    coor  = fl_table{:, 1:2};
    names = fl_table{:, 3};

end

