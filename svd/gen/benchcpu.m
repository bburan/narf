
disp('Benchmarching cpu speed:');

disp('Small memory usage...');
tic
N=100;
for ii=1:N,
   r=rand(100);
   [u,s,v]=svd(r);
end
x=toc;
fprintf('100 x 100 svd: %0.4f sec\n',x./N);

disp('Large memory usage...');
tic
N=1;
for ii=1:N,
   r=rand(1000);
   [u,s,v]=svd(r);
end
x=toc;
fprintf('1000 x 1000 svd: %0.4f sec\n',x./N);

