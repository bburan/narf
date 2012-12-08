
REPS=15;
fn=tempname;
%fn=['/auto/k1/david',fn];

for ii=1:REPS,
   d=rand(1000,5000);
   tic
   fid=fopen(fn,'wb');
   fwrite(fid,d,'double');
   %save(fn,'d');
   fprintf('Write: %.4f\n',toc);
   
   tic
   fid=fopen(fn,'rb');
   d2=fread(fid,size(d),'double');
   %d2=load(fn);
   fprintf('Read: %.4f\n',toc);
   fprintf('Diff: %.5f\n',sum(abs(d(:)-d2(:))));
   delete([fn]);
   
end


