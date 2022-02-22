function [] = cellarray2plots(data)
%CELLARRAY2PLOTS plot the time series data
%
%   example: cellarray2plots(data)
%
%   Input: a cell array 
%          Nx1 cells -> struct.var  -> Mx1 cells -> struct.data (double)
%                       struct.name (group)         struct.name (glacier)
%                       struct.time (datetime array)
%
%   Output: a figure that plots all time series labeled with names

    N = numel(data);
    dim = factor(N);
    if numel(dim) == 1 % it is a prime number
        dim = [1, dim]; % add a 2nd dim of 1
    end

    % start to plot
    % First grid the figure, meaning we define the subplot dimensions
    N_col = 7; % 7 subplots per row; this is fixed
    N_row = 3; % we start by guessing 3 rows
    while N_col*N_row < N
        N_row = N_row + 1; % add another row
    end
    
    N_subplot = N_col*N_row;
    dim = [N_row, N_col];
    
    
    fig1 = figure('Position', [50, 50, 1500, 800]);
    for i = 1:N_subplot
        
        if i > N % we don't have that many data, so just make an empty plot and go
            subplot(dim(1), dim(2), i);
            axis off
            continue
        else
            subplot(dim(1), dim(2), i)
        end
        % number of glacier flow lines we queried 
        M = numel(data{i,1}.var);
        time = data{i,1}.time;
        model_name = data{i,1}.name;
        
        % new loop: iterate over each glacier
        for j = 1:M
            glacier_name = data{i,1}.var{j,1}.name;
            rng(j*10)
            plot(time, data{i,1}.var{j,1}.data,...
                'Color', rand(1,3))
            hold on
        end
        hold off
        % force x axis limit
        xlim([time(1) time(end)])
        %legend
        xlabel('Time')
        ylabel('Thickness')
        title(model_name)
        
    end
        
end

