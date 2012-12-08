function I=quick_mi(x,y,xn,yn);

if ~exist('xn','var'),
   xn=10;
end
if ~exist('yn','var'),
   yn=xn;
end

xbins=linspace(min(x)-eps,max(x),xn+1);
ybins=linspace(min(y)-eps,max(y),yn+1);

% count the number of each discretized x
nmiu=zeros(xn,1);
for ii=1:xn,
   nmiu(ii)=sum(x>xbins(ii) & x<=xbins(ii+1));
end

% discretize y
yd=zeros(size(y));
for ii=1:yn,
   yd(find(y>ybins(ii) & y<=ybins(ii+1)))=ii-1;
end

% for each discete value of x, find all the yd
spk=zeros(1,1,max(nmiu),xn);
for ii=1:xn,
   spk(1,1,1:nmiu(ii),ii)=yd(find(x>xbins(ii) & x<=xbins(ii+1)));
end


bias=1; %quadratic extrapolation

% remove any bin that has less than three samples
fkeep=find(nmiu>3);
nmiu=nmiu(fkeep);
spk=spk(:,:,:,fkeep);

%mutual information. Plugin estimation
I=hr(spk,nmiu,bias)-hrs(spk,nmiu,bias);
