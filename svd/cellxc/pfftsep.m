function sepidx=pfftsep(H);

%H=H.*(H>0);

tsfIR=pfft2sf(H,'pfft');
Xmax=size(tsfIR,1);
obincount=24;
sfbincount=Xmax/2;

[tsfgr,obins,sfbins]=sf2gr(tsfIR,obincount,sfbincount);

mor=mean(tsfgr,2);

[u,s,v]=svd(tsfgr);
if u(:,1)'*mor >0,
   or0=u(:,1).*s(1);
   sf0=v(:,1);
else
   or0=-u(:,1).*s(1);
   sf0=-v(:,1);
end

tsfgr0=or0*sf0';

tsfIR1=gr2sf(tsfgr);
tsfIR0=gr2sf(tsfgr0);

sepidx=1-sqrt(sum((tsfIR1(:)-tsfIR0(:)).^2))./...
       sqrt(sum((tsfIR1(:)).^2));

return

mm=max(abs(tsfIR(:)));
mm0=max(abs(tsfIR1(:)));



figure(1)
clf
subplot(2,3,1);
imagesc(tsfIR,[-mm mm]);
axis image
subplot(2,3,2);
imagesc(tsfIR1,[-mm0 mm0]);
axis image
subplot(2,3,3);
imagesc(tsfIR0,[-mm0 mm0]);
axis image

subplot(2,3,6);
imagesc(abs(tsfIR0-tsfIR1),[-mm0 mm0]);
axis image

subplot(2,3,4);
imagesc(flipud(tsfgr'));
axis image
subplot(2,3,5);
imagesc(flipud(tsfgr0'));
axis image


