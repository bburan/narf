% kvafitnl.m
%
% specialized support code for various stages of kernvsatt.m
%
% created SVD 1/14/04 - extracted from xcfitatt.m 
%
% must be defined: USESIGMOIDS, attcount, tstrf
%                  linpred, respforfit (same length)
% affects tstrf
% creates nlinpred
%

nlcount=size(tstrf,1);
nlinpred=zeros(size(linpred,1),nlcount);

% set up stuff for fitting (relevant pred/obs time bins
% and associated att codes)
tokidx=find(sum(~isnan(respforfit(:,2:end)),2));
epred=linpred(tokidx,1);
fitresp=respforfit(tokidx,1);

if OUTNLMODE==2,
   % use hinge4 as output nl
   
   % nlidx=1 :  H same, dc/gain/thr same
   % nlidx=2 :  H same, dc float / gain/thr same
   % nlidx=3 :  H same, gain float / dc/thr same
   % nlidx=2 :  H same, thr float / dc/gain same
   % nlidx=5 :  H same, dc/gain/thr float
   % nlidx=6 :  H diff, dcg diff (dcg already removed for each att state)   
   
   % find baseline (all att) sigmoid parms
   tfitparms=fithinge4(epred,fitresp,0);
   for fitidx=1:attcount-1
      tstrf(1,fitidx).nltype='hinge4';
      tstrf(1,fitidx).nlparms=tfitparms;
   end
   nlinpred(tokidx,1)=hinge4(tfitparms,epred);
   
   spcount=3;
   for spidx=1:spcount,
      filler=ones(size(tokidx,1),spcount);
      for fitidx=1:attcount-1,
         filler(find(~isnan(respforfit(tokidx,fitidx+1))),spidx)=fitidx;
      end
      
      tfitparms=fithinge4(epred,fitresp,0,filler);
      
      for fitidx=1:attcount-1,
         useidx=[1:spidx-1 spidx+fitidx-1 ...
                 (spidx+attcount-1):(spcount+attcount-2)]';
         
         tstrf(spidx+1,fitidx).nltype='hinge4';
         tstrf(spidx+1,fitidx).nlparms=tfitparms(useidx);
         
         atokidx=find(~isnan(respforfit(tokidx,fitidx+1)));
         nlinpred(atokidx,spidx+1)=hinge4(tfitparms(useidx),...
                                           epred(atokidx));
      end
   end
   
   % nlidx=5: single kernel, float full output nl with attention
   for fitidx=1:attcount-1,
      attidx=fitidx+1;
      atokidx=find(~isnan(respforfit(:,attidx)));
      afitresp=respforfit(atokidx,attidx);
      epred=linpred(atokidx,1);
      fitparms=fithinge4(epred,afitresp);
      nlinpred(atokidx,5)=hinge4(fitparms,epred);
      tstrf(5,fitidx).nltype='hinge4';
      tstrf(5,fitidx).nlparms=fitparms;
   end
   
   % nlidx=6: output nl for separate att kernels. float dc/g
   % fit for each attidx
   for fitidx=1:attcount-1,
      attidx=fitidx+1;
      atokidx=find(~isnan(respforfit(:,attidx)));
      afitresp=respforfit(atokidx,attidx);
      epred=linpred(atokidx,attidx);
      fitparms=fithinge4(epred,afitresp);
      nlinpred(atokidx,6)=hinge4(fitparms,epred);
      tstrf(6,fitidx).nltype='hinge4';
      tstrf(6,fitidx).nlparms=fitparms;
   end
   
elseif OUTNLMODE==1,
   % ie, OUTNLMODE=1... use sigmoid NL
   
   % nlidx=1: fit dc and gain to all att states
   [sigparms]=fitsigmoid(epred,fitresp);
   nlinpred(tokidx,1)=sigmoid(sigparms,epred);
   for fitidx=1:attcount-1,
      tstrf(1,fitidx).nltype='sigmoid';
      tstrf(1,fitidx).nlparms=sigparms;
   end
   
   for nlidx=2:5,
      forcevalues=tstrf(1,1).nlparms;
      forcevalues(nlidx-1)=nan;
      
      for fitidx=1:attcount-1,
         attidx=fitidx+1;
         atokidx=find(~isnan(respforfit(:,attidx)));
         afitresp=respforfit(atokidx,attidx);
         epred=linpred(atokidx,attidx);
         
         sigparms=fitsigmoid(epred,afitresp,0,forcevalues);
         nlinpred(atokidx,nlidx)=sigmoid(sigparms,epred);
         tstrf(nlidx,fitidx).nltype='sigmoid';
         tstrf(nlidx,fitidx).nlparms=sigparms;
      end
   end
   
   % nlidx=6: output nl for separate att kernels. float dc/g
   % fit for each attidx
   for fitidx=1:attcount-1,
      attidx=fitidx+1;
      atokidx=find(~isnan(respforfit(:,attidx)));
      afitresp=respforfit(atokidx,attidx);
      epred=linpred(atokidx,attidx);
      %plot(linpred(atokidx,1),epred,'.');
      dcgparms=fitdcgain(epred,afitresp);
      nlinpred(atokidx,6)=dcgain(dcgparms,epred);
      tstrf(6,fitidx).nltype='dcgain';
      tstrf(6,fitidx).nlparms=dcgparms;
   end
elseif OUTNLMODE==1,
   % ie, OUTNLMODE=0... use dcgain
   
   % nlidx=1 :  H same, dc & gain same
   % nlidx=2 :  H same, dc float / gain same
   % nlidx=3 :  H same, dc same / gain float
   % nlidx=4 :  H same, dc & gain float
   % nlidx=5 :  H diff, single dcg (dcg already removed for each att state)
   % nlidx=6 :  H diff, dcg diff (dcg already removed for each att state)   
   attcode=zeros(size(tokidx));
   filler=ones(size(tokidx));
   for fitidx=1:attcount-1,
      attcode(find(~isnan(respforfit(tokidx,fitidx+1))))=fitidx;
   end
   
   % nlidx=1: fit dc and gain to all att states
   [dcgparms,beta0]=fitdcgain(epred,fitresp);
   nlinpred(tokidx,1)=dcgain(dcgparms,epred);
   for fitidx=1:attcount-1,
      tstrf(1,fitidx).nltype='dcgain';
      tstrf(1,fitidx).nlparms=dcgparms;
   end
   
   % nlidx=2:
   dcgparms=fitdcgain(epred,fitresp,[attcode filler]);
   nlinpred(tokidx,2)=dcgain(dcgparms,epred,[attcode filler+attcount-1]);
   for fitidx=1:attcount-1,
      tstrf(2,fitidx).nltype='dcgain';
      tstrf(2,fitidx).nlparms=[dcgparms(fitidx); dcgparms(attcount)];
   end
   
   % nlidx=3:
   dcgparms=fitdcgain(epred,fitresp,[filler attcode]);
   nlinpred(tokidx,3)=dcgain(dcgparms,epred,[filler attcode+1]);
   for fitidx=1:attcount-1,
      tstrf(3,fitidx).nltype='dcgain';
      tstrf(3,fitidx).nlparms=[dcgparms(1); dcgparms(fitidx+1)];
   end
   
   % nlidx=4:
   dcgparms=fitdcgain(epred,fitresp,[attcode attcode]);
   nlinpred(tokidx,4)=dcgain(dcgparms,epred,...
                             [attcode attcode+attcount-1]);
   for fitidx=1:attcount-1,
      tstrf(4,fitidx).nltype='dcgain';
      tstrf(4,fitidx).nlparms=...
          [dcgparms(fitidx); dcgparms(fitidx+attcount-1)];
   end
   
   % nlidx=5: single kernel, float full output nl with attention
   for fitidx=1:attcount-1,
      attidx=fitidx+1;
      atokidx=find(~isnan(respforfit(:,attidx)));
      afitresp=respforfit(atokidx,attidx);
      epred=linpred(atokidx,1);
      fitparms=fithinge4(epred,afitresp);
      nlinpred(atokidx,5)=hinge4(fitparms,epred);
      tstrf(5,fitidx).nltype='hinge4';
      tstrf(5,fitidx).nlparms=fitparms;
   end
   
   if 0, % old local att
   % nlidx=5: output nl for separate att kernels. single dc/g fit
   afitresp=respforfit(tokidx,1);
   epred=zeros(size(afitresp));
   for attidx=2:attcount,
      atokidx=find(~isnan(respforfit(tokidx,attidx)));
      epred(atokidx)=linpred(tokidx(atokidx),attidx);
   end
   
   dcgparms=fitdcgain(epred,afitresp);
   nlinpred(tokidx,5)=dcgain(dcgparms,epred);
   
   for fitidx=1:attcount-1,
      tstrf(5,fitidx).nltype='dcgain';
      tstrf(5,fitidx).nlparms=dcgparms;
   end
   end
   
   % nlidx=6: output nl for separate att kernels. float dc/g
   % fit for each attidx
   for fitidx=1:attcount-1,
      attidx=fitidx+1;
      atokidx=find(~isnan(respforfit(:,attidx)));
      afitresp=respforfit(atokidx,attidx);
      epred=linpred(atokidx,attidx);
      %plot(linpred(atokidx,1),epred,'.');
      dcgparms=fitdcgain(epred,afitresp);
      nlinpred(atokidx,6)=dcgain(dcgparms,epred);
      tstrf(6,fitidx).nltype='dcgain';
      tstrf(6,fitidx).nlparms=dcgparms;
   end
elseif OUTNLMODE==3,
   % ie, OUTNLMODE=3... fair dcgain
   
   % nlidx=1 :  H same, dc & gain same
   % nlidx=2 :  H same, dc float / gain same
   % nlidx=3 :  H same, dc & gain float
   % nlidx=4 :  H & gain rand, dc att
   % nlidx=4 :  H rand, gain att (dc fixed)
   % nlidx=6 :  H att, (dcg fixed)   
   attcode=zeros(size(tokidx));
   filler=ones(size(tokidx));
   for fitidx=1:attcount-1,
      attcode(find(~isnan(respforfit(tokidx,fitidx+1))))=fitidx;
   end
   
   % nlidx=1: fit dc and gain to all att states
   [dcgparms,beta0]=fitdcgain(epred,fitresp);
   nlinpred(tokidx,1)=dcgain(dcgparms,epred);
   for fitidx=1:attcount-1,
      tstrf(1,fitidx).nltype='dcgain';
      tstrf(1,fitidx).nlparms=dcgparms;
   end
   
   % nlidx=2:  global H, dc att
   dcgparms=fitdcgain(epred,fitresp,[attcode filler]);
   nlinpred(tokidx,2)=dcgain(dcgparms,epred,[attcode filler+attcount-1]);
   for fitidx=1:attcount-1,
      tstrf(2,fitidx).nltype='dcgain';
      tstrf(2,fitidx).nlparms=[dcgparms(fitidx); dcgparms(attcount)];
   end
   
   % nlidx=3:  global H, dcg att
   dcgparms=fitdcgain(epred,fitresp,[attcode attcode]);
   nlinpred(tokidx,3)=dcgain(dcgparms,epred,...
                             [attcode attcode+attcount-1]);
   for fitidx=1:attcount-1,
      tstrf(3,fitidx).nltype='dcgain';
      tstrf(3,fitidx).nlparms=...
          [dcgparms(fitidx); dcgparms(fitidx+attcount-1)];
   end
   
   % nlidx=4:  att H, global dcg
   epred=zeros(size(linpred(:,1)));
   for fitidx=1:attcount-1,
      attidx=fitidx+1;
      atokidx=find(~isnan(respforfit(:,attidx)));
      epred(atokidx)=linpred(atokidx,attidx);
   end
   epred=epred(tokidx);
   
   [dcgparms,beta0]=fitdcgain(epred,fitresp);
   nlinpred(tokidx,4)=dcgain(dcgparms,epred);
   for fitidx=1:attcount-1,
      tstrf(4,fitidx).nltype='dcgain';
      tstrf(4,fitidx).nlparms=dcgparms;
   end
   
   % nlidx=5:  att Hd, global g
   dcgparms=fitdcgain(epred,fitresp,[attcode filler]);
   nlinpred(tokidx,5)=dcgain(dcgparms,epred,[attcode filler+attcount-1]);
   for fitidx=1:attcount-1,
      tstrf(5,fitidx).nltype='dcgain';
      tstrf(5,fitidx).nlparms=[dcgparms(fitidx); dcgparms(attcount)];
   end
   
   % nlidx=6: output nl for separate att kernels. float dc/g
   % fit for each attidx
   for fitidx=1:attcount-1,
      attidx=fitidx+1;
      atokidx=find(~isnan(respforfit(:,attidx)));
      afitresp=respforfit(atokidx,attidx);
      epred=linpred(atokidx,attidx);
      dcgparms=fitdcgain(epred,afitresp);
      
      %nlinpred(atokidx,4)=dcgain(dcgparms,epred);
      %nlinpred(atokidx,5)=dcgain(dcgparms,epred);
      nlinpred(atokidx,6)=dcgain(dcgparms,epred);
      
      %tstrf(4,fitidx).nltype='dcgain';
      %tstrf(4,fitidx).nlparms=dcgparms;
      
      %tstrf(5,fitidx).nltype='dcgain';
      %tstrf(5,fitidx).nlparms=dcgparms;
      
      tstrf(6,fitidx).nltype='dcgain';
      tstrf(6,fitidx).nlparms=dcgparms;
   end
end



