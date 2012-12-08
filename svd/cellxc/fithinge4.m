% function [fitparms,beta0]=fithinge4(stim,resp,attcode);
%
% dcgparms=[dc g]
% y= dc + (g * x)
%
% beta0 = initial guess at dcgparms
% attcode allows certain parameters to float with attention state.
%
% NOTE: assumes attcode ranges from 1...max for each column
%
function [fitparms,beta0]=fithinge4(stim,resp,showfit,attcode);

if ~exist('showfit','var'),
   showfit=0;
end

if ~exist('attcode','var'),
   attcode=ones(length(stim),3);
   attcounts=[1 1 1];
else
   % assume attcode ranges from 1...max for each column
   attcounts=max(attcode);
end

% shift reference to align with position in beta vector
attcode(:,2)=attcode(:,2)+attcounts(1);
attcode(:,3)=attcode(:,3)+attcounts(1)+attcounts(2);

r0=resp-mean(resp);
r1=stim-mean(stim);
d1=sum(r1.^2);

if d1==0,
   % no modulation in predicted response
   fitparms=[ones(attcounts(1),1) .* mean(resp); 
             zeros(attcounts(2),1);
             zeros(attcounts(3),1)];
   beta0=fitparms;
   
   return
end

scf=sum(r0.*r1)./d1;

beta0=[ones(attcounts(1),1) .* min(resp);
       ones(attcounts(2),1)  .* scf; 
       ones(attcounts(3),1)  .* mean(stim)./5  ];
lb=[ones(attcounts(1),1)  .* -inf;
    ones(attcounts(2),1)  .* 0.0;
    ones(attcounts(3),1)  .* -inf];
ub=[ones(attcounts(1),1)  .* inf;
    ones(attcounts(2),1)  .* inf;
    ones(attcounts(3),1)  .* inf   ];

fitopt=optimset('Display','off');

fitparms=lsqcurvefit('hinge4',beta0,stim,resp,lb,ub,fitopt,attcode);


%[beta0 fitparms]


% plot fit over data
if showfit,
   hold off
   ptshow=min([length(stim) 500]);
   
   scatter(stim(1:ptshow),resp(1:ptshow),'.');
   xx=linspace(min(stim),max(stim),50);
   hold on
   plot(xx,hinge4(beta0,xx),'r--');
   plot(xx,hinge4(fitparms,xx));
   hold off
   drawnow
   title(sprintf('hinge fit: %.2f %.2f %.2f',fitparms));
   legend('init','fit');
end
