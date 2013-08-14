

rawdata=mysql(['SELECT gDataRaw.*,gPenetration.numchans',...
               ' FROM gDataRaw,gCellMaster,gPenetration',...
               ' WHERE gDataRaw.masterid=gCellMaster.id',...
               ' AND gCellMaster.penid=gPenetration.id',...
               ' AND gCellMaster.siteid >= "lim042"',...
               ' AND gCellMaster.siteid <= "lim047"',...
               ' AND gDataRaw.runclass="TOR"',...
               ' AND not(gDataRaw.bad)']);
%               ' AND gCellMaster.siteid >= "lim042"',...
%               ' AND gCellMaster.siteid <= "lim043"',...
%               ' AND gCellMaster.siteid >= "lim036"',...
%               ' AND gCellMaster.siteid <= "lim020"',...
%               ' AND gCellMaster.siteid >= "dnb034"',...
%               ' AND gCellMaster.siteid <= "dnb036"',...

%chrange=[3 4 5 7 8 10 11 12 13 17 18 19 23 24 28 30];
%chrange=1:32;
chrange=11;
%chrange=[9 10 11 12 20 25];
%chrange=1:9;
fnstr='elec';
%fnstr='chan';

for ii=1:length(rawdata),
    parmfile=[rawdata(ii).resppath rawdata(ii).parmfile];
    ff=figure;
    rr=ceil(sqrt(length(chrange)));
    cc=ceil(length(chrange)./rr);
    for jj=1:length(chrange),
        ch=chrange(jj);
        tfile=[rawdata(ii).resppath 'tmp/' ...
               strrep(rawdata(ii).parmfile,'.m','') ...
               sprintf('.001.1.%s%d.sig4.mat',fnstr,ch)];
        tfile=strrep(tfile,'/auto/data/Danube','/auto/data/daq/Danube');
        if exist(tfile,'file'),
            load(tfile);
            
            sfigure(ff);
            subplot(rr,cc,jj);
            %spikematrix=spikematrix-repmat(mean(spikematrix),...
            %                              [size(spikematrix,1) 1]);
            if size(spikematrix,2)>200,
                plot(spikematrix(:,round(linspace(1,size(spikematrix,2),200))));
            else
                plot(spikematrix);
            end
            axis tight
            aa=axis;
            hold on
            plot(aa(1:2),[0 0],'k--');
            hold off
        end
        
        if jj==1,
            ht=title(basename(parmfile));
        else
            ht=title(sprintf('%d (%d spk)',ch,size(spikematrix,2)));
        end
        set(ht,'Interpreter','none');
        axis off
        drawnow
    end
    
end
