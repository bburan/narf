function [resp,cdata,fdata]=xcgetbootset(params,bootidx,times,resp0,stim,raster);

if ~exist('raster','var'),
   raster=[];
end

   filelens=times(1).stop-times(1).start+1;
   filecount=max(find(filelens>1));
   filelens=filelens(1:filecount);

   totlen=sum(filelens);
   clen=[0; cumsum(filelens)];
   
   bootcount=params.fitboot;
   bootstep=filelens./bootcount;
   spacer=diff(params.maxlag);
      % segregate out the exploratory, fit and validation data
      % in order to keep things balanced, take the same fraction
      % from each file included in the stim and resp matrices
      fdata.stim=[];
      fdata.resp=[];

      % if predfrac>0 use the same validation data for each
      % bootstrap. this is done in the standard section far
      % below. otherwise create anew validation set during each boot
      % iteration
      cdata.stim=[];
      cdata.resp=[];
      cdata.raster=[];
      cdata.fromidx=[];
      cdata.toidx=[];

      resp=resp0;
      for fidx=1:filecount,

         % pull out bootfrac of the time bins for each fitting and
         % validation
         if params.predfrac==0,
            valstart=round(clen(fidx)+(bootidx-1)*bootstep(fidx)+1);
            valstop=round(clen(fidx)+bootidx*bootstep(fidx));
         else
            valstart=1;
            valstop=0;
         end

         bidx2=bootidx+round(bootcount/2);
         if params.fitfrac==0,

            % still set these values to exclude them from analysis
            % for proper jackknifing
            fitstart=1;
            fitstop=0;
            %fitstart=round(clen(fidx)+mod(bidx2-1,bootcount)*bootstep(fidx)+1);
            %fitstop=round(clen(fidx)+(mod(bidx2-1,bootcount)+1)*...
            %              bootstep(fidx));
         else

            % change this to deal with a variable fit range
            %fitrange
            fitstart=round(clen(fidx)+mod(bidx2-1,bootcount)*bootstep(fidx)+1);
            fitstop=round(clen(fidx)+(mod(bidx2-1,bootcount)+1)*...
                          bootstep(fidx));
         end

         fprintf('BOOT %2d/%2d: fidx=%d fit: %6d-%6d  val: %6d-%6d / %6d\n',...
                 bootidx,bootcount,fidx,fitstart,fitstop,...
                 valstart,valstop,totlen);

         % reset resp and nan out the ranges reserved for fit and pred
         resp(fitstart:fitstop,:)=nan;
         resp(valstart:valstop,:)=nan;
         
         if isfield(params,'respfiltercmd') && ...
               strcmp(params.respfiltercmd,'resptakefracs'),
            respfracs=params.respfilterparms{1};
            resprange=(clen(fidx)+1):clen(fidx+1);
            
            goodidx=find(~isnan(resp(resprange,1)))+clen(fidx);
            goodlen=length(goodidx);
            shift=fitstop;

            resp(resprange,2:end)=nan;
            for ff=2:length(respfracs),
               keepidx=1:round(length(goodidx).*respfracs(ff));
               keepidx=mod(keepidx+shift-1,goodlen)+1;
               %keepidx=round(linspace(1,length(goodidx),...
               %                       length(goodidx).*respfracs(ff)));
               resp(goodidx(keepidx),ff)=resp(goodidx(keepidx),1);
            end
         end

         if params.fitfrac>0,
            s0=min([params.maxlag(2) fitstart-clen(fidx)-1]);
            e0=min([-params.maxlag(1) clen(fidx+1)-fitstop]);

            fdata.stim=cat(1,fdata.stim,stim(fitstart-s0:fitstop+e0,:));

            tr=resp0(fitstart-s0:fitstop+e0,:);
            tr([1:s0 length(tr)-e0+1:length(tr)],:)=nan;
            fdata.resp=cat(1,fdata.resp,tr);
         end

         if params.predfrac==0,
            s0=min([params.maxlag(2) valstart-clen(fidx)-1]);
            e0=min([-params.maxlag(1) clen(fidx+1)-valstop]);

            cdata.fromidx=[cdata.fromidx valstart:valstop];
            cdata.toidx=[cdata.toidx ...
                         ((s0+1):(valstop-valstart+s0+1))+length(cdata.resp)];
            
            cdata.stim=cat(1,cdata.stim,stim(valstart-s0:valstop+e0,:));
            
            tr=resp0(valstart-s0:valstop+e0,:);
            tr([1:s0 length(tr)-e0+1:length(tr)],:)=nan;
            cdata.resp=cat(1,cdata.resp,tr);
            
            if ~isempty(raster),
               traster=raster(valstart-s0:valstop+e0,:);
               traster([1:s0 length(tr)-e0+1:length(tr)],:)=nan;
               cdata.raster=cat(1,cdata.raster,traster);
            end
            
         end
      end

 
