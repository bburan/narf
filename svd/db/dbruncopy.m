% function [runidx,fileidx]=dbruncopy(runidxin,sresfile,newbatch,totlens);
%
% copy an entire run entry (tRunData, tRunFile) with a new resfile
% specified in sresfile
%
function [runidx,fileidx]=dbruncopy(runidxin,sresfile,newbatch,totlens);

dbopen;

sql=['SELECT * FROM tRunData where id=',num2str(runidxin)];
rundata=mysql(sql);

sql=['SELECT * FROM tRunFile where rundataid=',num2str(runidxin),...
     ' ORDER BY usecode desc'];
filedata=mysql(sql);

if exist('newbatch'),
   if length(newbatch)>1,
      disp('newbatch must be a scalar!');
      return
   end
end
if ~exist('totlens'),
   totlens=0;
end

for lenidx=1:length(totlens),
   f=fieldnames(rundata);
   snames='INSERT INTO tRunData (';
   svalues=') VALUES (';
   bfirst=1;
   
   for ii=1:length(f)
      if not(strcmp(f{ii},'id')),
         t=getfield(rundata,f{ii});
         if strcmp(f{ii},'resfile'),
            if length(totlens)>0 & sum(totlens)>0,
               if totlens(lenidx)==0,
                  sv=['"',sresfile,'.all','"'];
               else
                  sv=['"',sresfile,'.',sprintf('%.5d',(totlens(lenidx))),'"'];
               end
            else
               sv=['"',sresfile,'"'];
            end
         elseif isempty(t),
            sv='null';
         elseif isnumeric(t),
            sv=num2str(t);
         else
            sv=['"',t,'"'];
         end
         
         if bfirst,
            snames=[snames,f{ii}];
            svalues=[svalues,sv];
            bfirst=0;
         else
            snames=[snames,',',f{ii}];
            svalues=[svalues,',',sv];
         end
      end
   end

   sql=[snames,svalues,');']
   [res,aff,runidx]=mysql(sql);
   %runidx=lenidx;

   % if newbatch is specified, modify it:
   if exist('newbatch'),
      sql=['UPDATE tRunData set complete=0,batch=',num2str(newbatch),...
           ' WHERE id=',num2str(runidx)];
      mysql(sql);
   end
   
   fileidx=zeros(length(filedata),1);
   nn=1;
   totlensofar=0;
   while nn<=length(filedata) & ...
             (totlensofar<totlens(lenidx) | totlens(lenidx)==0),
      f=fieldnames(filedata);
      snames='INSERT INTO tRunFile (';
      svalues=') VALUES (';
      bfirst=1;
      
      for ii=1:length(f)
         if strcmp(f{ii},'respstop') & filedata(nn).usecode==0,
            resplen=filedata(nn).respstop-filedata(nn).respstart+1;
            
            if totlens(lenidx)>0 & totlensofar+resplen>totlens(lenidx),
               resplen=totlens(lenidx)-totlensofar;
               respstop=filedata(nn).respstart+resplen-1;
            else
               respstop=filedata(nn).respstop;
            end
            totlensofar=totlensofar+resplen;
            
            sv=num2str(respstop);
            
            if bfirst,
               snames=[snames,f{ii}];
               svalues=[svalues,sv];
               bfirst=0;
            else
               snames=[snames,',',f{ii}];
               svalues=[svalues,',',sv];
            end
         
         elseif not(strcmp(f{ii},'id')),
            t=getfield(filedata(nn),f{ii});
            if strcmp(f{ii},'rundataid'),
               sv=num2str(runidx);
            elseif isempty(t),
               sv='null';
            elseif isnumeric(t),
               sv=num2str(t);
            else
               sv=['"',t,'"'];
            end
            
            if bfirst,
               snames=[snames,f{ii}];
               svalues=[svalues,sv];
               bfirst=0;
            else
               snames=[snames,',',f{ii}];
               svalues=[svalues,',',sv];
            end
         end
      end
      
      sql=[snames,svalues,');'];
      [res,aff,fileidx(nn)]=mysql(sql);
      %fileidx(nn)=nn;
      
      nn=nn+1;
   end
end

