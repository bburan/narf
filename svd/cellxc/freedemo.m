% freedemo.m
%
% demo of some useful functions for visualizing and loading v1/v4 data


% list the free view reverse correlation runs
%dblist(27);  % list v1 review runs
%dblist(30);  % list v4 freeview runs


% loading v1 demo

% figure out runidx for particular cellid and batch
cellid='R166C';
batch=27;
sql=['SELECT id FROM tRunData',...
     ' WHERE cellid="',cellid,'" AND batch=',num2str(batch)];
rundata=mysql(sql);
runidx=rundata(1).id;

decorrspace=3;
predtype=0;    % 0 for V1, 2 for V4
bstep=1;
fitidx=3;

%movmatchres(runidx,decorrspace,predtype);  

[sfIR,sfIRsum,sfIRall,kernfile]=loadkern(runidx,decorrspace,...
                                         predtype,bstep,fitidx);

ttime=squeeze(sum(sum(sum(sfIR,1),2),4));
tt=min(find(ttime==max(ttime)));

msfIR=squeeze(sfIR(:,:,tt,:));

figure(1);
clf
plot((0:(length(ttime)-1))*14,ttime);

movdecorrres(msfIR,0,[cellid,' peak resp'],2,4);



% loading v4 demo

% figure out runidx for particular cellid and batch
cellid='m0048';
batch=30;  % freview batch
sql=['SELECT id FROM tRunData',...
     ' WHERE cellid="',cellid,'" AND batch=',num2str(batch)];
rundata=mysql(sql);
runidx=rundata(1).id;

decorrspace=3;
predtype=2;    % 0 for V1, 2 for V4
bstep=0;       % don't use step
fitidx=3;

%movmatchres(runidx,decorrspace,predtype);  

[sfIR,sfIRsum,sfIRall,kernfile]=loadkern(runidx,decorrspace,...
                                         predtype,bstep,fitidx);
% don't do this
%ttime=squeeze(sum(sum(sum(sfIR,1),2),4));
%tt=min(find(ttime==max(ttime)));
% instead, always use +25 through +75 bin
tt=6;

msfIR=squeeze(sfIR(:,:,tt,:));

figure(1);
clf
plot((0:(length(ttime)-1))*14,ttime);

movdecorrres(msfIR,0,[cellid,' peak resp'],2,4);
