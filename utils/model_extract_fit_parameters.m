function [res,filecodes]=model_extract_fit_parameters(modelpath);
    
    global STACK XXX;
    
    load_model(char(modelpath));
    res={};
    filecodes=XXX{1}.filecodes;
    for ii=1:length(STACK),
        if isfield(STACK{ii}{1},'fit_fields'),
            fitfields=STACK{ii}{1}.fit_fields;
            for jj=1:length(fitfields),
                for kk=1:length(STACK{ii}),
                    res{ii,kk}.(fitfields{jj})= ...
                        STACK{ii}{kk}.(fitfields{jj});
                end
            end
        end
    end
 