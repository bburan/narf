% function [rcomp]=compact_raster_matrix3(rraw)
%
% takes a sinlge trial raster matrix with gaps (eg, from only showing
% segments of a movie on individual trials but preseriving the
% individual trial results) and collapses it so that there's still
% only a single trial per column
%
% input: rraw - raster matrix (time X trial)... each entry contains
% the number of spikes occuring in that time bin during that
% trial. invalid bins should be filled with "nan"
%
% output: rcomp - compacted raster (time X effective trials)
%
% created SVD 2003 - modified from BV's compact_raster_matrix2.m
% 
function [rcomp]=compact_raster_matrix3(rraw)

rok=rraw>=0;
rnz=sum(rok,2);
rcount=min(rnz(find(rnz>0)));

rcomp=ones(size(rraw,1),rcount).*nan;

for rr=1:size(rraw,2);
   rok=find(rraw(:,rr)>=0);
   rrr=sum(isnan(rcomp(rok,:)),1);
   insidx=min(find(rrr));
   
   if ~isempty(insidx),
      rcomp(rok,insidx)=rraw(rok,rr);
   end
end

