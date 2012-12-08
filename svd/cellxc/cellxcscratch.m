
clear xcparms

RECALC=1;

global BATQUEUEID
BATQUEUEID=[];
cellid='r0284';
speed=60;
runclassid=2;
fmt='pixel';
respfmtcode=0;

%'speed',speed,
cellfiledata=dbgetscellfile('cellid',cellid,...
                            'speed',speed,...
                            'respfmtcode',respfmtcode,...
                            'runclassid',runclassid,'fmt',fmt);
sql=['SELECT * FROM gCellMaster WHERE id=',...
     num2str(cellfiledata(1).masterid)];
celldata=mysql(sql);

xcparms.stimfiles={};
xcparms.respfiles={};
filecount=0;
optidx=[];
for ii=1:length(cellfiledata),
   if findstr('-opt-',cellfiledata(ii).stimfile),
      optidx=[optidx ii];
      fprintf('%d Found opt-file: %s\n',ii,cellfiledata(ii).stimfile);
      
   elseif findstr('-synth-',cellfiledata(ii).stimfile),
      optidx=[optidx ii];
      fprintf('%d Found synth-file: %s\n',ii,cellfiledata(ii).stimfile);
      
   elseif findstr('conf',cellfiledata(ii).stimfile),
      optidx=[optidx ii];
      fprintf('%d Found conf-file: %s\n',ii,cellfiledata(ii).stimfile);
      
   else
      
      fprintf('%d Found exp-file: %s\n',ii,cellfiledata(ii).stimfile);
      filecount=filecount+1;
      xcparms.stimfiles{filecount}=[cellfiledata(ii).stimpath,...
                    cellfiledata(ii).stimfile];
      xcparms.respfiles{filecount}=[cellfiledata(ii).path,...
                    cellfiledata(ii).respfile];
   end
end

for ii=1:length(optidx),
   filecount=filecount+1;
   xcparms.stimfiles{filecount}=[cellfiledata(optidx(ii)).stimpath,...
                    cellfiledata(optidx(ii)).stimfile];
   xcparms.respfiles{filecount}=[cellfiledata(optidx(ii)).path,...
                    cellfiledata(optidx(ii)).respfile];
end


[ll,iconside]=imfileinfo(xcparms.stimfiles{1});
scalepix=round(iconside(1)./celldata.rfsize*8);

xcparms.cellid=cellid;

xcparms.stimloadcmd='loadimfile';
%xcparms.stimloadparms={12,17,12};

if 0,
   xcparms.stimloadparms={scalepix,8,16};
   xcparms.stimfiltercmd='movphasesep';
   xcparms.stimfilterparms={0,0,0,1,0};
   xcparms.kernfmt='fft';
elseif 1,
   xcparms.stimloadparms={scalepix,0,16};
   xcparms.stimfiltercmd='movpower';
   xcparms.stimfilterparms={0,0,1,1,0};
   xcparms.kernfmt='pfft';
elseif 0
   addpath /auto/k1/jess/ssa/
   xcparms.stimloadparms={scalepix,scalepix+1,16};
   xcparms.stimfiltercmd='sca';
   xcparms.stimfilterparms={};
   xcparms.kernfmt='space';
else
   xcparms.stimloadparms={scalepix,scalepix+1,16};
   xcparms.stimfiltercmd='';
   xcparms.stimfilterparms={};
   xcparms.kernfmt='space';
end

% for fv data

%xcparms.resploadcmd='resploadatt';
%xcparms.resploadparms={'psth',1,0,1};

xcparms.resampfmt=1;
xcparms.resampcount=20;
xcparms.outfile=['/auto/k5/david/tmp/' cellid 'out.mat'];
xcparms.sfsstep=6;
xcparms.sfscount=30;
xcparms.sffiltsigma=5;
xcparms.smoothtime=0;
xcparms.fitfrac=0.05;
xcparms.predfrac=0.05;
xcparms.fitboot=0;
xcparms.docellfit2=0;
xcparms.repexclude=0;
xcparms.nloutparm=4;

if RECALC,
   cellxcnodb(xcparms);
else
   xcresult(xcparms.outfile);
end

disp('stopped before optanal');
keyboard

if length(optidx)>0,
   % optinat analysis!
   cellfiledata=dbgetscellfile('cellid',cellid,...
                               'respfmtcode',respfmtcode,...
                               'runclassid',runclassid,'fmt',fmt);
   
   load(xcparms.outfile);
   
   nlidx=2;
   
   ir=sum(strf(nlidx).h .* (strf(nlidx).h>0), 1);
   ir=ir(1:8);
   maxlag=min(find(max(ir)==ir))-1;
   
   for tpredfile=optidx,
      respfile=[cellfiledata(tpredfile).path,...
                cellfiledata(tpredfile).respfile];
      tr=feval(params.resploadcmd,respfile,params.resploadparms{:});
      stimfile=[cellfiledata(tpredfile).stimpath,...
                cellfiledata(tpredfile).stimfile];
      %mov=loadimfile(stimfile);
      smov=loadimfile(stimfile,0,0,64,65,32,1);
      
      cparms=params;
      cparms.stimfiles={stimfile};
      cparms.respfiles={respfile};
      cparms.stimcrfs=4;
      tpredfile=1;
      tpredstartframe=1;
      tpredstopframe=size(tr,1);
      [cdata.stim,cdata.resp]=xcloadstimresp(tpredfile,tpredstartframe,...
                                             tpredstopframe,cparms);
      optres=xcval(strf(nlidx),params,cdata);
      
      pp=optres.mod_psth{1}(1:tpredstopframe);
      rr=optres.act_resp{1}(1:tpredstopframe);
      pp(find(isnan(pp)))=0;
      rr(find(isnan(rr)))=0;
      
      %keyboard
      
      hir=fft(ir);
      hir=hir./abs(hir);
      fir=real(ifft(hir));
      decp=conv(pp,flipud(fir'));
      decr=conv(rr,flipud(fir'));
      
      %decp=decp(1:(end-length(fir)+1));
      %decr=decr(1:(end-length(fir)+1));
      decp=decp(length(fir):end);
      decr=decr(length(fir):end);
      
      decp=decp./mean(decp).*mean(pp);
      decr=decr./mean(decr).*mean(rr);
      
      [junk,ppi]=sort(-decp);
      [junk,rri]=sort(-decr);
      
      
      
      %
      % plot psth and sorted responses
      %
      figure(1);
      clf
      
      subplot(3,1,1);
      plot(rr,'k--');
      hold on
      plot(pp,'r');
      hold off
      legend('obs','pred');
      title(sprintf('%s optinat analysis. predcorr=%.2f',...
                    params.cellid,optres.predxc));
      
      subplot(3,1,2);
      plot(decr,'k--');
      hold on
      plot(decp,'r');
      hold off
      legend('deconv obs','deconv pred');
      
      subplot(3,1,3);
      plot(decp(ppi),'k--');
      hold on
      plot(decr(ppi),'r');
      hold off
      legend('pred ranked','obs');
      
      %plot(decr(length(ir):end));
      %hold on
      %%plot(rr,'r--');
      %hold off
      %legend('dec obs','obs');
      
      %showkern(log(cat(3,cdata.stim(ppi(1:10),:)',...
      %                 cdata.stim(rri(1:10),:)')),params.kernfmt);
      
      %
      % compute angle stuff
      %
      tbincount=maxlag;
      h=sum(strf(nlidx).h(:,1:maxlag+1),2);
      spacecount=length(h);
      mS=strf(nlidx).mS;
      tstim=cdata.stim'-repmat(mS,1,size(cdata.stim,1));
      %tstim=[zeros(spacecount,tbincount) tstim(:,1:(end-tbincount))];
      
      %h=h./std(tstim,1,2);
      %tstim=tstim ./ repmat((std(tstim,1,2)),1,size(tstim,2));
      
      linpred=tstim'*h;
      
      tstim(find(abs(h)==0),:)=0;
      
      snorm=sqrt(sum(tstim.^2,1))';
      snorm(find(snorm==0))=1;
      angle=acos(linpred./norm(h(:))./snorm) *180/pi;
      
      figure(2);
      clf
      
      scatter(angle,decp-decr);
      xlabel('angle');
      ylabel('pred-obs resp');
      
      if 0,
         plot3(linpred,angle,decr,'o');
         xlabel('linpred');
         ylabel('stim . h angle');
         zlabel('obs resp');
         grid on
         
         tp=sqrt(decp+min(decp));
         tr=sqrt(decr+min(decr));
         
         scatter(tp,angle,decr*20+10);
         xlabel('linpred');
         ylabel('stim . h angle');
         
         scatter(tp,tr,sqrt((angle-min(angle)+1))*20);
         xlabel('sqrt(linpred)');
         ylabel('sqrt(act)');
      end
      
      %
      % matched array stuff.
      %
      prank=ones(size(ppi));
      prank(ppi)=(length(ppi):-1:1);
      rrank=ones(size(ppi));
      rrank(rri)=(length(rri):-1:1);
      
      decb=rrank; % (prank+rrank)./2;
      dece=(rrank-prank)./2;
      %decb=(decp+decr)./2;
      %dece=(decp-decr)./(decb+1);
      [junk,bbi]=sort(-decb);
      [junk,eei]=sort(-dece);
      
      DOPREDERR=0;
      DORANK=1;    % 0-plot response strength
                   % 1-plot ranks 
                   % 2-plot some weird error thing vs actual resp rank
      DOTRANS=0;   % ie, do pfft of images for display
      
      if DOTRANS,
         psize=16;
      else
         psize=32;
      end
      bpsize=30*psize;
      pplen=length(ppi);
      bigplot=uint8(ones(bpsize).*255);
      showcount=pplen;
      
      if DOPREDERR,
         dece=decr-decp;
         rmin=min([dece; decr]);
         tp=(dece-rmin);
         tr=(decr-rmin);
         trmin=0;
      else
         rmin=min([decp; decr]);
         tp=(decp-rmin);
         tr=(decr-rmin);
         trmin=0;
      end
      
      rrange=max([tp; tr])-rmin;
      
      figure(4);
      clf
      
      kernfmt=params.kernfmt;
      spacebincount=size(cdata.stim,2);
      for ii=showcount:-1:1,
         
         if DORANK==1,
            pidx=ppi(ii);
            xpos=round((bpsize-psize)./pplen*(ii-1));
            ypos=round((bpsize-psize)./pplen*(find(rri==ppi(ii))-1));
            rstr='rank';
         elseif DORANK==2,
            pidx=bbi(ii);
            %xpos=round((bpsize-psize) ./ rrange .* (tp(pidx)-trmin));
            %ypos=round((bpsize-psize) ./ rrange .* (tr(pidx)-trmin));
            
            xpos=round((bpsize-psize)./pplen*(ii-1));
            ypos=round((bpsize-psize)./pplen*(find(eei==bbi(ii))-1));
            rstr='rank2';
         else
            pidx=ppi(ii);
            xpos=round((bpsize-psize) ./ rrange .* (tp(pidx)-trmin));
            ypos=round((bpsize-psize) ./ rrange .* (tr(pidx)-trmin));
            rstr='resp';
         end
         
         if DOTRANS & (strcmp(kernfmt,'fft') | strcmp(kernfmt,'pfft')),
            
            fim=cdata.stim(pidx,:);
            %fim=cdata.stim(pidx,:)'-mS;
            %fim=fim-h./norm(h).^2 .* decp(pidx).^2;
            
            if strcmp(kernfmt,'fft'),
               phasecount=4;
            else
               phasecount=1;
            end
            
            chancount=spacebincount/phasecount;
            Xmax=sqrt(chancount*2);
            [cfilt,cfiltconj,cmask,cmaskconj]=gencfilt(Xmax,Xmax);
            tsf=zeros(Xmax,Xmax);
            tsf(cfilt)=sum(reshape(fim,chancount,phasecount),2);
            tsf(cfiltconj)=tsf(cfilt);
            tsf=reshape(tsf,Xmax,Xmax);
            psize=size(tsf,1);
            im=tsf;
            %im=sqrt(tsf-min(tsf(:)));
            im(find(im(:)<0))=0;
            im=im.^(0.33);
            if max(im(:))>0,
               im=im./max(im(:)).*512;
            end
         else
            im=smov(:,:,pidx);
         end
         
         bigplot((xpos+1):(xpos+psize),(ypos+1):(ypos+psize))=im;
      end
      if DORANK,
         hc=imagesc(linspace(1,showcount,bpsize),...
                    linspace(1,showcount,bpsize),bigplot);
      else
         hc=imagesc(linspace(rmin,rmin+rrange,bpsize),...
                    linspace(rmin,rmin+rrange,bpsize),bigplot);
      end
      ha=get(hc,'Parent');
      tl=get(ha,'YTickLabel');
      
      xlabel(['actual ' rstr]);
      ylabel(['pred ' rstr]);
      axis image
      colormap(gray);
      title(sprintf('cell %s, optinat tested with %s kernel',...
                    params.cellid,params.kernfmt));
      set(gcf,'PaperPosition',[0.25 0.25 8 10.5],...
              'PaperOrientation','portrait');
      
      if 0,
         figure(3);
         clf
         showcount=20;
         for ii=1:20,
         
         subplot(4,showcount/2,ii);
         imagesc(smov(:,:,ppi(ii)),[0 255]);
         axis image
         axis off
         title(sprintf('pred: %d',ppi(ii)));
         
         subplot(4,showcount/2,ii+showcount);
         imagesc(smov(:,:,rri(ii)),[0 255]);
         axis image
         axis off
         title(sprintf('act: %d',rri(ii)));
         end
      colormap(gray)
      end
      
      if length(optidx)>1,
         disp('pausing before next optfile');
         keyboard;
      end
   end  
   
end