


cc=repmat(round(linspace(1,64,6)),[2 1]);
cc=repmat(cc(:)',[16 1 1 16]);
jj=jet;

for ii=1:16,
   shiftset=mod(ii+(0:3),16)+1;
   cc(shiftset,4:9,:,17-ii)=cc(shiftset,9:-1:4,:,17-ii);
end
imwrite(cc,jj,'/auto/p1/svd/html/music/images/lrainbow_an.gif',...
        'LoopCount',inf,'DelayTime',0.1);
imwrite(flipdim(cc,2),jj,'/auto/p1/svd/html/music/images/rrainbow_an.gif',...
        'LoopCount',inf,'DelayTime',0.1);
imwrite(permute(cc,[2 1 3 4]),jj,...
        '/auto/p1/svd/html/music/images/rurainbow_an.gif',...
        'LoopCount',inf,'DelayTime',0.1);

imwrite(permute(flipdim(cc,1),[2 1 3 4]),jj,...
        '/auto/p1/svd/html/music/images/lurainbow_an.gif',...
        'LoopCount',inf,'DelayTime',0.1);


for ii=1:size(cc,4),
   imagesc(cc(:,:,:,ii));
   pause(0.1);
end
