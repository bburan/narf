% function emailres(mailto,scommand,subject)
%
% created SVD 2/01
%
% execute matlab command(s) contained in mailcommand and email the
% output to mailto (via unix mail command)
%
% mailto - (default 'stephen@socrates.berkeley.edu' if empty)
%
function emailres(mailto,scommand,subject)

DEFMAILTO='svdavid@berkeley.edu';

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

unix(['/usr/sbin/sendmail ',mailto,' < ',fn]);
delete(fn);


