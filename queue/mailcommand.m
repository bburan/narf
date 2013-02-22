% function mailcommand(mailto,scommand,subject)
%
% created SVD 2/01
%
% execute matlab command(s) contained in mailcommand and email the
% output to mailto (via unix mail command)
%
% mailto - (default 'stephen@socrates.berkeley.edu' if empty)
%
function mailcommand(mailto,scommand,subject)

DEFMAILTO='stephen@socrates.berkeley.edu';

if isempty(mailto),
   mailto=DEFMAILTO;
end

fn=tempname;

diary(fn);
disp(['To: ',mailto]);
if ~exist('subject','var'),
   subject=['scommand: ',scommand];
end
disp(['Subject: ',subject]);
eval(scommand);
diary off

%unix(['mail ',mailto,' < ',fn]);
%if checkbic==0,
   disp('sending mail via jlg');
   unix(['/usr/sbin/sendmail -C/auto/k1/david/.sendmail.cf ',...
         mailto,' < ',fn]);
%else
%   disp('sending mail via bic/millennium');
%   unix(['mail -s "',subject,'" ',mailto,' < ',fn]);
%end
delete(fn);

