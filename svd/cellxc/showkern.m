%function showkern(H,kernfmt,iconside,titles,contours,binsize,frcount);
%
% 
function showkern(H,kernfmt,iconside,titles,contours,binsize,frcount);

global GCOLORMAP

if ~exist('kernfmt','var'),
   kernfmt='pfft+4';
end
if ~exist('titles','var'),
   titles={};
end
if ~exist('contours','var'),
   contours=0;
end
if ~exist('binsize','var'),
   binsize=14;
end
LINEPLOT=0;

spacebincount=size(H,1);
tbincount=size(H,2);
kcount=size(H,3);
if strcmp(kernfmt,'fft'),
   kernfmt='pfft+4';
end

if strcmp(kernfmt,'fft') | strcmp(kernfmt,'sfft') | ...
      strcmp(kernfmt,'pfft') | ...
      ~isempty(findstr(kernfmt,'pfft-')) | ...
      ~isempty(findstr(kernfmt,'pfft+')) | ...
      strcmp(kernfmt,'lin2'),
   
   if strcmp(kernfmt,'fft') | strcmp(kernfmt,'lin2'),
      phasecount=4;
   elseif length(kernfmt)>4 & ...
         (kernfmt(end-1)=='+' | kernfmt(end-1)=='-'),
      phasecount=str2num(kernfmt(end));
   else
      phasecount=1;
   end
   
   chancount=spacebincount/phasecount;
   Xmax=sqrt(chancount*2);
   
   if Xmax~=round(Xmax),
      chancount=(spacebincount-1)/phasecount;
      Xmax=sqrt(chancount*2);
      H=H(1:spacebincount-1,:,:,:,:);
   end
   
   [cfilt,cfiltconj,cmask,cmaskconj]=gencfilt(Xmax,Xmax);
   
   if kernfmt(end-1)=='-',
      tsf=zeros(Xmax*Xmax,tbincount,kcount);
      tsf(cfilt,:,:)=squeeze(sum(reshape(H,chancount,phasecount,...
                                         tbincount,kcount),2));
      tsf(cfiltconj,:,:)=tsf(cfilt,:,:);
      tsf=reshape(tsf,Xmax,Xmax,tbincount,kcount);
      phasecount=1;
   else
      % don't average over phase
      phx=ceil(sqrt(phasecount));
      phy=ceil(phasecount/phx);
      
      tsf=zeros(Xmax*Xmax,tbincount,kcount,phx*phy);
      for ii=1:phasecount,
         tsf(cfilt,:,:,ii)=H((1:chancount)+(ii-1)*chancount,:,:);
      end
      tsf(cfiltconj,:,:,:)=tsf(cfilt,:,:,:);
      
      tsf([1 2 Xmax],:)=repmat(tsf(3,:),[3 1]);
      
      tsf=reshape(tsf,Xmax,Xmax,tbincount,kcount,phx,phy);
      
      tsf=reshape(permute(tsf,[1 5 2 6 3 4]),...
                  Xmax*phx,Xmax*phy,tbincount,kcount);
   end
   
elseif strcmp(kernfmt,'fftgr') | strcmp(kernfmt,'pfftgr'),
   
   if strcmp(kernfmt,'fftgr'),
      phasecount=4;
   else
      phasecount=1;
   end
   chancount=spacebincount/phasecount;
   Xmax=sqrt(chancount*2);
   
   obins=linspace(0,180,15+1);
   obins=obins(1:end-1);
   sfbins=linspace(1,Xmax/2,Xmax/2);
   if strcmp(kernfmt,'fftgr'),
      tsf=sf2gr(H,length(obins),length(sfbins),0,0,'pfft+4');
   else
      tsf=sf2gr(H,length(obins),length(sfbins),0,0,'pfft');
   end
   
   tsf=permute(tsf,[2 1 3 4]);
   tsf=flipdim(tsf,1);
   
   if 0
   obincount=8;
   sfbincount=8;
   if strcmp(kernfmt,'fftgr'),
      phasecount=4;
      tsf=squeeze(fft2gr(H,obincount,sfbincount,phasecount));
      
      tsf=flipdim(tsf,2);
      tsf=reshape(permute(tsf,[2 4 1 3 5]),sfbincount*phasecount,...
                  obincount,tbincount,kcount);
   else
      phasecount=1;
      tsf=squeeze(fft2gr(H,obincount,sfbincount,phasecount));
      tsf=flipdim(tsf,2);
      tsf=permute(tsf,[2 1 3 4]);
   end
   end
elseif strcmp(kernfmt,'logfft2'),
   phasecount=4;
   chancount=spacebincount/phasecount;
   Xmax=sqrt(chancount);
   tsf=reshape(H,Xmax,Xmax,phasecount,tbincount,kcount);
   
   if kcount==1,
      tsf=permute(tsf,[1 2 4 3]);
   else
      tsf=squeeze(sum(tsf,3));
   end
   
elseif strcmp(kernfmt,'wav') | strcmp(kernfmt,'wg'),
   %reshape to x/or(fine) X y/sf(fine) X time X kernel
   tsf=reshape(H,[iconside tbincount kcount]);
   tsf=permute(tsf,[4 1 3 2 5 6]);
   tsf=reshape(tsf,iconside(4)*iconside(1),iconside(3)*iconside(2),...
               tbincount,kcount);

elseif strcmp(kernfmt,'winfft'),
   Xmax=iconside(1);
   Ymax=iconside(2);
   specbincount=iconside(3);
   H=reshape(H,[Xmax Ymax specbincount tbincount kcount]);
   
   sfcount=sqrt(specbincount*2);
   [cfilt,cfiltconj]=gencfilt(sfcount,sfcount);
   tsf=zeros(sfcount,Xmax,sfcount,Ymax,tbincount,kcount);
   tsf0=zeros(sfcount*2,tbincount,kcount);
   for ix=1:Xmax,
      for iy=1:Ymax,
         tsf0(cfilt,:,:)=reshape(H(ix,iy,:,:,:),specbincount,tbincount,kcount);
         tsf0(cfiltconj,:,:)=reshape(H(ix,iy,:,:,:),specbincount,...
                                     tbincount,kcount);
         tsf(:,ix,:,iy,:,:)=reshape(tsf0,sfcount,1,sfcount,1,...
                                    tbincount,kcount);
      end
   end
   tsf=reshape(tsf,sfcount*Xmax,sfcount*Ymax,tbincount,kcount);
   iconside=[Xmax Ymax sfcount sfcount];
elseif strcmp(kernfmt,'mfilt'),
   tsf=reshape(H,[iconside tbincount kcount]);
   tsf=permute(tsf,[1 3 2 4 5]);
   tsf=reshape(tsf,iconside(1)*iconside(3),iconside(2),...
               tbincount,kcount);
   
elseif strcmp(kernfmt,'orsf'),
   if length(iconside)==1,
      iconside=[iconside 1];
   end
   tsf=reshape(H,[iconside tbincount kcount]);   
   
elseif strcmp(kernfmt,'gr'),
   tsf=reshape(H,[iconside tbincount kcount]);
   
elseif strcmp(kernfmt,'fourier'),
   phasecount=2;
   chancount=spacebincount/phasecount;
   Xmax=sqrt(chancount*2);
   [cfilt,cfiltconj,cmask,cmaskconj]=gencfilt(Xmax,Xmax);
   
   if 0, % old,m complex space
      tsf=zeros(Xmax*Xmax,tbincount,kcount);
      tsf(cfilt,:,:)=abs(reshape(H,chancount,tbincount,kcount));
      tsf(cfiltconj,:,:)=tsf(cfilt,:,:);
      tsf=reshape(tsf,Xmax,Xmax,tbincount,kcount);
   else
      tsf=zeros(Xmax*Xmax*phasecount,tbincount,kcount);
      tsf(cfilt,:,:)=H(1:chancount,:,:);
      tsf(cfiltconj,:,:)=tsf(cfilt,:,:);
      tsf(cfilt+spacebincount,:,:)=H((1:chancount)+chancount,:,:);
      tsf(cfiltconj+spacebincount,:,:)=tsf(cfilt+spacebincount,:,:);
      
      tsf=reshape(tsf,Xmax,Xmax,2,tbincount,kcount);
      tsf=permute(tsf,[1 3 2 4 5]);
      tsf=reshape(tsf,Xmax*2,Xmax,tbincount,kcount);
   end
   
elseif strcmp(kernfmt,'lin2'),
   phasecount=4;
   chancount=spacebincount/phasecount;
   Xmax=sqrt(chancount*2);
   [cfilt,cfiltconj,cmask,cmaskconj]=gencfilt(Xmax,Xmax);
   
   tsf=zeros(Xmax*Xmax*2,2,tbincount,kcount);
   tsf(cfilt,1,:,:)=H(1:chancount,:,:);
   tsf(cfiltconj,1,:,:)=H(1:chancount,:,:);
   tsf(cfilt+chancount*2,1,:,:)=H((1:chancount)+chancount,:,:);
   tsf(cfiltconj+chancount*2,1,:,:)=H((1:chancount)+chancount,:,:);
   tsf(cfilt,2,:,:)=H((1:chancount)+chancount*2,:,:);
   tsf(cfiltconj,2,:,:)=H((1:chancount)+chancount*2,:,:);
   tsf(cfilt+chancount*2,2,:,:)=H((1:chancount)+chancount*3,:,:);
   tsf(cfiltconj+chancount*2,2,:,:)=H((1:chancount)+chancount*3,:,:);
   %keyboard
   tsf=reshape(tsf,Xmax,Xmax,2,2,tbincount,kcount);
   tsf=permute(tsf,[1 3 2 4 5 6]);
   tsf=reshape(tsf,Xmax*2,Xmax*2,tbincount,kcount);
   
elseif strcmp(kernfmt,'spect'),
   % plot like standard spectrogram plots
   iconside=spacebincount;
   tsf=reshape(H,iconside,tbincount,1,kcount);
   tsf=permute(tsf,[1 2 3 4]);
   tsf=flipdim(tsf,1);
   tbincount=1;
   if size(tsf,1)==1,
      LINEPLOT=1;
   elseif size(tsf,1)>3,
      tsf=imresize(tsf,2,'bilinear');
   elseif size(tsf,1)>1,
      tsf=tsf;
   else
      tsf=imresize(tsf,2,'bilinear');
      tsf=permute(tsf,[1 2 4 3]);
   end
   %smooth=[100 250];
   %tsf = interpft(interpft(tsf,smooth(2),2),smooth(1),1);
elseif strcmp(kernfmt,'space') | strcmp(kernfmt,'pix') | ...
      strcmp(kernfmt,'pixel') | strcmp(kernfmt,'sca'),
   if ~exist('iconside','var') | prod(iconside)~=size(H,1),
      Xmax=round(sqrt(spacebincount));
      if Xmax*Xmax==size(H,1),
         tsf=reshape(H,Xmax,Xmax,tbincount,kcount);
      else
         tsf=zeros(Xmax*Xmax,tbincount,kcount);
         tsf(1:size(H,1),:,:)=H;
         tsf=reshape(tsf,Xmax,Xmax,tbincount,kcount);
      end
   elseif length(iconside)==1,
      Xmax=round(sqrt(spacebincount));
      iconside=[Xmax spacebincount./Xmax];
      tsf=reshape(H,[iconside tbincount kcount]);
   else
      tsf=reshape(H,[iconside tbincount kcount]);
   end
end

% actually do the display

if ~exist('frcount','var'),
   frcount=15;
end
if frcount>tbincount,
   frcount=tbincount;
end
if not(exist('orientation','var')),
   orientation=1;
end

%clf reset
sfscount=size(tsf,4);
for r=1:sfscount,
   a=(tsf(:,:,:,r));
   astd=std(a(:));
   amax=max(abs(a(:)));
   amin=-amax;
   
   if isnan(H(1,1,r)),
      tfrcount=0;
   else
      tfrcount=frcount;
   end
   
   for fr=1:tfrcount,
      if orientation==1,
         h=subplot(sfscount,frcount,fr+frcount*(r-1));
      else
         h=subplot(frcount,sfscount,(fr-1)*sfscount+r);
      end
      
      if LINEPLOT,
          plot(a(:,:,fr));
      elseif amin~=amax,
          imagesc(a(:,:,fr),[amin,amax]);
      else
          imagesc(zeros(size(a,1),size(a,2)));
      end
      
      if strcmp(kernfmt,'fftgr') | strcmp(kernfmt,'pfftgr') | strcmp(kernfmt,'space'),
         xx=-2:2;
         gsf=exp(-(xx./1.0).^2/2);
         gsf=gsf./sum(gsf(:));
         gor=exp(-(xx./1.2).^2/2);
         gor=(gor./sum(gor(:)))';
         
         if contours,
            ta=a(:,:,fr);
            ta=conv2(ta,gsf,'same');
            ta=cconv2(ta,gor);
            tas=size(ta);
            
            slev1=1.0;
            slev2=2.0;
            hold on
            contour(ta,[-astd*slev2  -astd*slev2],'k-');
            contour(ta,[-astd*slev1  -astd*slev1],'k--');
            contour(ta,[ astd*slev1   astd*slev1],'k--');
            contour(ta,[ astd*slev2   astd*slev2],'k-');
            hold off
         end
      elseif strcmp(kernfmt,'wav') | strcmp(kernfmt,'wg') | ...
            strcmp(kernfmt,'winfft'),
         hold on
         for ii=1:(iconside(1)-1),
            hl=plot([0.5 size(a,2)+0.5],[ii*iconside(4) ii*iconside(4)]+0.5,...
                    'k','LineWidth',0.25);
         end
         for ii=1:(iconside(2)-1),
            hl=plot([ii*iconside(3) ii*iconside(3)]+0.5,[0.5 size(a,1)+0.5],...
                    'k','LineWidth',0.25);
         end
         
         if 0,
            % obins go verticle to horizontal to vertica..
            %keyboard
            obincount=iconside(3);
            sfbincount=iconside(4);
            obins=linspace(0,pi,obincount+1);
            obins=obins(1:end-1);
            obins=exp(i.*obins(:).*2);
            for xx=1:iconside(1),
               for yy=1:iconside(2),
                  cf=a(:,:,fr);
                  patch=cf(((xx-1)*obincount)+1:xx*obincount,...
                          ((yy-1)*sfbincount)+1:yy*sfbincount);
                  [uu,ss,vv]=svd(patch);
                  oc=uu(:,1);
                  sc=vv(:,1);
                  A=ss(1,1);
                  if sum(oc)<0,
                     oc=-oc;
                     sc=-sc;
                  end
                  om=angle(mean(oc.*obins))/2;
                  os=min(find(sc==max(sc)));
                  
                  %kathleen, look here
                  cx=mean([((xx-1)*obincount)+1 xx*obincount]);
                  cy=mean([((yy-1)*sfbincount)+1 yy*sfbincount]);
                  dx=(obincount/2)*cos(om);
                  dy=-(obincount/2)*sin(om);
                  hl=plot([cy-dy cy+dy],[cx-dx cx+dx],'k-');
                  set(hl,'LineWidth',os,...
                         'Color',[1 1 1].*(1-max(abs(patch(:)))./amax));
               end
            end
         end
         
         
         hold off
      elseif ~isempty(findstr(kernfmt,'pfft+')),
         phasecount=str2num(kernfmt(end));
         phx=ceil(sqrt(phasecount));
         phy=ceil(phasecount/phx);
         
         axx=axis;
         hold on
         for ii=1:(phx-1),
            hl=plot([0.5 size(a,2)+0.5],[ii*Xmax ii*Xmax]+0.5,...
                    'k','LineWidth',0.25);
         end
         for ii=1:(phy-1),
            hl=plot([ii*Xmax ii*Xmax]+0.5,[0.5 size(a,2)+0.5],...
                    'k','LineWidth',0.25);
         end
         hold off
         axis(axx);
      else
         
         if 0 & contours,
            
            xx=-2:2;
            gsf=exp(-(xx./0.75).^2/2);
            gsf=gsf./sum(gsf(:));
            gor=exp(-(xx./0.75).^2/2);
            gor=(gor./sum(gor(:)))';
            
            ta=conv2(a(:,:,fr),gsf,'same');
            ta=conv2(ta,gor,'same');
            tas=size(ta);
            
            if strcmp(kernfmt,'fft') | strcmp(kernfmt,'pfft')
               xc=ceil((tas(1)+1)/2);
               yc=ceil((tas(2)+1)/2);
               ta(xc,yc)=0;
            end
               
            hold on
            contour(ta,[-astd*3.0  -astd*3.0],'k-');
            contour(ta,[-astd*1.5  -astd*1.5],'k--');
            contour(ta,[ astd*1.5   astd*1.5],'k--');
            contour(ta,[ astd*3.0   astd*3.0],'k-');
            hold off
         end
      end
      
      if size(a,1)>1 & size(a,2)>1 & ~strcmp(kernfmt,'fftgr') & ~strcmp(kernfmt,'spect'),
         axis image
      end
      
      set(h,'YTickLabel',[]);
      set(h,'XTickLabel',[]);
      if r==sfscount & orientation==1,
         if fr==frcount,
            xlabel('Latency (ms)','FontSize',11);
         else
            xlabel(sprintf('%d',round((fr-1)*binsize+binsize/2)),...
                   'FontSize',11);
         end
      end
      if r==1 & fr==1 & (strcmp(kernfmt,'wav') | strcmp(kernfmt,'wg')),
         xlabel('x/or','FontSize',11);
         ylabel('y/sf','FontSize',11);
      end
      if fr==round(frcount/2) & r<=length(titles),
         ht=title(sprintf('%s',titles{r}));
         set(ht,'FontSize',11,'Interpreter','none');
      end
   end
end

if ~isempty(GCOLORMAP),
   colormap(GCOLORMAP);
elseif strcmp(kernfmt,'space') | strcmp(kernfmt,'pix'),
   colormap(gray);
elseif strcmp(kernfmt,'spect'),
   %colormap(redblue);
   colormap('default');
else
   colormap(redblue);
   %colormap(blueyellow);
end

%spkcolormap=[0 0 1; 0.5 0.5 1; 1 1 1; 1 0.5 0.5; 1 0 0];
%spkcolormap=[0 0 1; 1 1 1; 1 1 1; 1 1 1; 1 0 0];
%colormap(spkcolormap);

% force it to fill entire landscape page
%colormap(gray);
%if orientation==2,
   set(gcf,'PaperOrientation','landscape','PaperPosition',[0.25 0.25 10.5 8]);
%else
%   set(gcf,'PaperOrientation','portrait','PaperPosition',[0.25 0.25 8 10.5]);
%end

drawnow;

