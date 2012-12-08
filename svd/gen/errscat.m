% function errscat(d1,d2,e1,e2);
%
%
function errscat(d1,d2,e1,e2);

hold on
for ii=1:length(d1),
   if ~isnan(d1(ii)) & ~isnan(d2(ii)),
      plot([d1(ii) d1(ii)],[d2(ii)-e2(ii) d2(ii)+e2(ii)]);
      plot([d1(ii)-e1(ii) d1(ii)+e1(ii)],[d2(ii) d2(ii)]);
   end
end
hold off
