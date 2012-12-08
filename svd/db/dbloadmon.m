% dbloadmon.m : monitor rc queue activity


reloadtime=5;
disp('Starting dbloadmon.m. Press ctrl-c to end.');

OVERSENSITIVE=1; % kill jobs on bic machines if load goes over some threshold
USERCOLOR=1;     % plot job progress for all users with entries in queue
SHOWALLJLG=1;

% to get load off remote computer: ssh -1 <host> cat /proc/loadavg
% (on ketamine :  ssh1 <host> cat /proc/loadavg)

while 1,
   
   [r,compnames]=dbgetload([],1,15,1);
   r(:,2)=-r(:,2);
   r(:,3)=-(r(:,3)>0);
   [r,ix]=sortrows(r,[3 2]);
   r(:,2)=-r(:,2);
   compnames={compnames{ix}};
   if SHOWALLJLG,
      activecomps=find(r(:,3) | r(:,4)==0 | r(:,4)==3);
   else
      activecomps=find(r(:,3));
   end
   if length(activecomps)==0,
      fprintf('\nNo active jobs. Exiting dbloadmon.\n');
      break
   end
   
   if OVERSENSITIVE,
      hiactive=find(r(:,4)==1 & r(:,2)>1.5 & r(:,3)>0);
      for ii=1:length(hiactive),
         fprintf('%s: Overload! Killing jobs on %s.\n',...
                 datestr(now),compnames{hiactive(ii)});
         sql=sprintf(['UPDATE tQueue set killnow=1',...
                      ' WHERE complete=-1 AND machinename="%s"'],...
                     [compnames{hiactive(ii)},'.bic']);
         mysql(sql);
      end
   end
   
   if USERCOLOR,
      sql=['SELECT count(id) as count,user,complete',...
           ' FROM tQueue GROUP BY user,complete'];
      userdata=mysql(sql);
      
      users={};
      completerange=[1 -1 2 0];
      complete=[];
      for ii=1:length(userdata),
         curcomplete=find(completerange==userdata(ii).complete);
         curuser=find(strcmp(users,userdata(ii).user));
         if isempty(curuser),
            users={users{:},userdata(ii).user};
            curuser=length(users);
         end
         complete(curuser,curcomplete)=userdata(ii).count;
      end
      usercount=length(users);
      compfrac=length(activecomps)./(length(activecomps)+length(users));
   end
   
   figure(1);
   clf
   colormap(winter);
   
   if USERCOLOR,
      subplot('position',[0.05 0.1+(1-compfrac)*0.85 0.9 compfrac*0.85]);
   end
   
   tr=r(activecomps,2);
   tjobs=-r(activecomps,3);
   rn=tr.*(tjobs==0);
   rg=tr.*(tr<1.1 & tjobs>0);
   rr=tr.*(tr>=1.1 & tjobs>0);
   
   barh([rn rg rr],'stacked');
   
   legend('inactive','low load','high load');
   
   for ii=1:length(activecomps),
      sname=sprintf('%s (%d)',compnames{activecomps(ii)},...
                    r(activecomps(ii),5));
      ht=text(r(activecomps(ii),2)+0.05,ii,sname,...
              'HorizontalAlignment','left','FontSize',16);
   end
   title('load per computer');
   
   if USERCOLOR,
      subplot('position',[0.05 0.05 0.9 (1-compfrac)*0.85]);
      barh(complete,'stacked');
      for ii=1:length(users),
         ht=text(sum(complete(ii,:))+0.05,ii,users{ii},...
                 'HorizontalAlignment','left','FontSize',16);
      end
      legend('complete','in prog','dead','not started',-1);
      title('users');
   end
   
   set(1,'WindowButtonDown','drawnow');
   pause(reloadtime);
   
end

