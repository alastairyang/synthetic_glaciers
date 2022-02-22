function keyword = combine_keys(varargin)
%COMINE_KEYS combine multiple keys and output a keyword argument for dir()

    N = nargin;
    keyword = [];
    for i = 1:N
        if isa(varargin{i}, 'string')
            varargin{i} = convertStringsToChars(varargin{i});
        end
        keyword = [keyword, '*', varargin{i}];
    end
    % crop out the beginning *
    keyword = convertCharsToStrings(keyword(2:end));
        
end

