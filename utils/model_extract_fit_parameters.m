function res=model_extract_fit_parameters(modelpath);
    
    global STACK XXX;
    
    load_model(char(modelpath));
    res={};
    for ii=1:length(STACK),
        if isfield(STACK{ii}{1},'fit_fields'),
            fitfields=STACK{ii}{1}.fit_fields;
            for jj=1:length(fitfields),
                res{ii}.(fitfields{jj})= ...
                    STACK{ii}{1}.(fitfields{jj});
            end
        end
    end
 