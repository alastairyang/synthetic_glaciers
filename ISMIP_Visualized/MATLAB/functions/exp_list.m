function list = exp_list(data)
%EXP_LIST Create a list of experiment names and their indices for future
%reference
%   input:
%       data: all experiment data in a cell array. Each cell is a structure
%   output:
%       list: a list in a table format. First column numerical index,
%       second column the strings

    N = numel(data);
    names = [];
    index = [];
    for i = 1:N
        if isa(data{i}.att.institution, 'char')
            name = convertCharsToStrings(data{i}.att.institution);
        else
            name = data.att.institution;
        end
        names = [names; name];
        index = [index; i];
    end
    
    % convert to a table
    list = table(index, names);
    
end

