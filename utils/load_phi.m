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

