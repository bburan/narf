% function [phi,phi_label,narfdata]=load_phi(cellid,batch,modelname);
%
% load all parameters that were fit for a model or set of models
% matching the (cellid,batch,modelname) search criteria in the
% NarfResults table.  
%
% inputs:
% cellid - baphy cellid (or search string)
% batch - baphy/narf data set number
% modelname - narf model name
%
% outputs:
% phi - a [parameter count X cell] matrix of fit parameter values,
%       where each column contains all parameters for one model fit
% phi_label - {parameter count X cell} cell array with a string
%             identifying the model and parameter number of each
%             value in phi.
% narfdata - [cell X 1] structure array with metadata for each model
%
% "%" is a wildcared, so this sequence of commands will load fit
% parameters for a bunch of cells fit with the same model: 
% >> modelname='env100_logfree_dep1pc_fir15_siglog100_fit05';
% >> batch=259;
% >> cellid='%';
% >> [phi,phi_label,narfdata]=load_phi(cellid,batch,modelname);
%
% created SVD 2014_05_23
%
function [phi,phi_label,narfdata]=load_phi(cellid,batch,modelname);

global STACK

sql=['SELECT * FROM NarfResults WHERE cellid like "',cellid,'"',...
     ' AND batch=',num2str(batch),...
     ' AND modelname like "',modelname,'"'];
narfdata=mysql(sql);

if isempty(narfdata),
    error('no (cellid,batch,modelname) match in NarfResults');
end

modelcount=length(narfdata);
phi_label={};
phi=[];

for midx=1:modelcount,
    tphi_label={};
    tphi=[];
    modelpath=char(narfdata(midx).modelpath);
    load_model(modelpath);
    for ii=1:length(STACK),
        if isfield(STACK{ii}{1},'fit_fields') && ...
                ~isempty(STACK{ii}{1}.fit_fields),
            for jj=1:length(STACK{ii}{1}.fit_fields),
                fn=STACK{ii}{1}.fit_fields{jj};
                phi1=STACK{ii}{1}.(fn)(:);
                tphi=cat(1,tphi,phi1);
                for kk=1:length(phi1),
                    tphi_label=cat(1,tphi_label,{...
                        sprintf('%d_%s_%s_%02d',ii,STACK{ii}{1}.name,fn,kk)});
                end
            end
        end
    end
    phi=cat(2,phi,tphi);
    phi_label=cat(2,phi_label,tphi_label);
end

