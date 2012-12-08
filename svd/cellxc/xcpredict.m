% function r=kernpredict(strf,mov)
%
% strf      - strf structure(s)
% mov       - time X spacecount : (linearized) movie matrix
%
% created SVD 11/19/04
%
function r=xcpredict(strf,mov)

strfcount=size(strf);
movlen=size(mov,1);
spacecount=size(mov,2);

r=zeros([movlen strfcount]);

for strfidx=1:prod(strfcount),
   tstim=mov'-repmat(strf(strfidx).mS,[1 movlen]);
   seplinpred=kernpredict(strf(strfidx).h,tstim,spacecount,0);
   linpred=sum(seplinpred,2);
   
   if strcmp(strf(strfidx).nltype,'none'),
      r(:,strfidx)=linpred;
      
   elseif strcmp(strf(strfidx).nltype,'post-sum thresh'),
      
      r(:,strfidx)=thresh(strf(strfidx).nlparms,linpred);
      
   elseif strcmp(strf(strfidx).nltype,'exp thresh'),
      
      r(:,strfidx)=expthresh(strf(strfidx).nlparms,linpred);
      
   elseif strcmp(strf(strfidx).nltype,'pre-sum thresh'),
      
      r(:,strfidx)=hinge2(strf(strfidx).nlparms,seplinpred);
   elseif strcmp(strf(strfidx).nltype,'suppnorm'),
      
      % figure out which channels are inhibitory and which are excitatory
      hp=find(strf(strfidx).hspace>0);
      hn=find(strf(strfidx).hspace<0);
      tr=[sum(seplinpred(:,hp),2) -sum(seplinpred(:,hn),2)];
      r(:,strfidx)=suppnorm(strf(strfidx).nlparms,tr);
   else
      % assume nltype is a string
      r(:,strfidx)=feval(strf(strfidx).nltype,...
                         strf(strfidx).nlparms,tmod_psth);
   end

end
