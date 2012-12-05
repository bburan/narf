function coefficients = stephen_fit();

mystim = permute(XXX{end}.dat.por022b02_p_TOR.ds_stim, [2 1 3]);
mystim = reshape(mystim, 30*600, 10);
myresp = XXX{end}.dat.por022b02_p_TOR.raw_respavg;
myresp = myresp';
myresp = myresp(:);

STRF = cellxcdataloaded(stim,resp,params);


coefficents = STRF.h ;