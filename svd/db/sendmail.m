% function sendmail(mailto,stext,subject)
%
% created SVD 2/01
%
% execute matlab command(s) contained in mailcommand and email the
% output to mailto (via unix mail command)
%
% mailto - (default 'stephen@socrates.berkeley.edu' if empty)
%
function sendmail(mailto,stext,subject)

DEFMAILTO='stephen@kamzik.org';

if isempty(mailto),
   mailto=DEFMAILTO;
end

fn=tempname;

diary(fn);
disp(['To: ',mailto]);
disp(['From: root@bhangra.isr.umd.edu']);
disp(['Subject: ',subject]);
disp(stext);
diary off

[w,s]=unix(['/usr/sbin/sendmail ',mailto,' < ',fn]);
delete(fn);


