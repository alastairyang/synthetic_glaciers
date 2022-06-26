%% Model Syn_0
index = '0';
whole_suite(index);

%% Model Syn_00
index = '00';
whole_suite(index);

%% Model Syn_01
index = '01';
whole_suite(index);

%% Function
function [] = whole_suite(index)

    % all vars
    runname = 'meltrates_rheoB_fric';
    model_out = run_models(index, 't','all');
    if strcmp('0',index)
        md = model_out.model_0_t;
        path = ['results/syn_',index,'/',runname,'.mat'];
        save(path, 'md')
    end
    if strcmp('00',index)
        md = model_out.model_00_t;
        path = ['results/syn_',index,'/',runname,'.mat'];
        save(path, 'md')
    end
    if strcmp('01',index)
        md = model_out.model_01_t;
        path = ['results/syn_',index,'/',runname,'.mat'];
        save(path, 'md')
    end
    
    
    % fric
    runname = 'fric';
    model_out = run_models(index, 't','fric');
    if strcmp('0',index)
        md = model_out.model_0_t;
        path = ['results/syn_',index,'/',runname,'.mat'];
        save(path, 'md')
    end
    if strcmp('00',index)
        md = model_out.model_00_t;
        path = ['results/syn_',index,'/',runname,'.mat'];
        save(path, 'md')
    end
    if strcmp('01',index)
        md = model_out.model_01_t;
        path = ['results/syn_',index,'/',runname,'.mat'];
        save(path, 'md')
    end
    
    % meltrates
    runname = 'meltrates';
    model_out = run_models(index, 't','meltrates');
    if strcmp('0',index)
        md = model_out.model_0_t;
        path = ['results/syn_',index,'/',runname,'.mat'];
        save(path, 'md')
    end
    if strcmp('00',index)
        md = model_out.model_00_t;
        path = ['results/syn_',index,'/',runname,'.mat'];
        save(path, 'md')
    end
    if strcmp('01',index)
        md = model_out.model_01_t;
        path = ['results/syn_',index,'/',runname,'.mat'];
        save(path, 'md')
    end

    % rheoB
    runname = 'rheoB';
    model_out = run_models(index, 't','rheoB');
    if strcmp('0',index)
        md = model_out.model_0_t;
        path = ['results/syn_',index,'/',runname,'.mat'];
        save(path, 'md')
    end
    if strcmp('00',index)
        md = model_out.model_00_t;
        path = ['results/syn_',index,'/',runname,'.mat'];
        save(path, 'md')
    end
    if strcmp('01',index)
        md = model_out.model_01_t;
        path = ['results/syn_',index,'/',runname,'.mat'];
        save(path, 'md')
    end

end