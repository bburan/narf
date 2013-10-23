function fitset=export_fit_parameters(batch, cellids, modelnames)

    function ret = getit(x, idx)
        if isempty(x)
            ret = nan;
        else
            ret = x(idx);
        end
    end
    
    function ret = getfitval(x)
        if ~isfield(x, 'r_fit') || ~isfield(x, 'r_test')
            ret = [nan nan];
        else
            ret = [x.r_fit x.r_test];
        end
    end

stack_extractor = @(x) x;
meta_extractor = @getfitval;
% stacks is modelnames X cellid
[stacks, metas, x0s,preds] = load_model_batch(batch, cellids, modelnames, ...
                                        stack_extractor, meta_extractor);
cellcount=length(cellids);
modelcount=length(modelnames);

fitset=struct();

for cc=1:cellcount,
    for mm=1:modelcount,
        fitset(cc,mm).cellid=cellids{cc};
        fitset(cc,mm).modelname=modelnames{mm};
        fitset(cc,mm).training_set=x0s{mm,cc}.training_set;
        fitset(cc,mm).test_set=x0s{mm,cc}.test_set;
        fitset(cc,mm).filecodes=x0s{mm,cc}.filecodes;
        
        fitset(cc,mm).r_test=preds{mm,cc}(1);
        fitset(cc,mm).r_fit=preds{mm,cc}(2);
        
        unique_codes = unique(x0s{mm,cc}.filecodes);
        
        for ii=1:length(stacks{mm,cc}),
            for jj=1:length(stacks{mm,cc}{ii}),
                if isfield(stacks{mm,cc}{ii}{jj},'fit_fields'),
                    ff=stacks{mm,cc}{ii}{jj}.fit_fields;
                    for ff=1:length(stacks{mm,cc}{ii}{jj}.fit_fields),
                        fname=stacks{mm,cc}{ii}{jj}.fit_fields{ff};
                        if length(stacks{mm,cc}{ii})>1,
                            stepname=sprintf('p%02d%s_%s_%s',ii,...
                                             stacks{mm,cc}{ii}{1}.name,...
                                             fname,unique_codes{jj});
                        else
                            stepname=sprintf('p%02d%s_%s',ii,...
                                             stacks{mm,cc}{ii}{1}.name,fname);
                        end
                        fitset(cc,mm).(stepname)=...
                            stacks{mm,cc}{ii}{jj}.(fname);
                    end
                end
            end
        end
        
   end
end

if nargout==0,
    [FileName,PathName,FilterIndex] = ...
        uiputfile(['narf_export_', num2str(batch), '.mat'],'Save to...');
    if ~FileName,
        disp('Canceled save.');
    else
        fprintf('Saving fitset to %s%s\n',PathName,FileName);
        save([PathName FileName],'fitset','cellids','modelnames');
    end
end


end

