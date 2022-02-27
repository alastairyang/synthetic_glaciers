function [] = sens_data_sampler(model_index, test_index, sens_data)
%SENS_DATA_SAMPLER This functions samples a subset of the sensitivity dataset 
% for simulations
%
% Input:
%       model_index[str]: model index
%       test_index[int] : test run index
%       sens_data[struc]: structure containing the data

    modelname = ['sens_', model_index];
    data = sens_data.(modelname);
    N_test = data.N_test;
    N_vars  = data.N_vars;

    [var_n, test_n] = ind2sub([N_vars, N_test], test_index);
    teststr = ['test_', num2str(test_n)];
    varnames = fieldnames(data);

    % if the variable data is a double array, meaning this forcing data is
    % not time-dependent (which comes in cell array form)
    % We turn this to a cell too
    if isequal(class(data.(varnames{var_n}).(teststr)), 'double')
        datacell = cell(1);
        datacell{1} = data.(varnames{var_n}).(teststr);
        sampled.data = datacell;
    else
        sampled.data = data.(varnames{var_n}).(teststr);
    end

    sampled.var = varnames{var_n};

    save('sampled vars/sampled.mat','sampled')
end

