function cfs = MakeErbCFs(mincf,maxcf,numchans)
% 
% MAKEERBCFS Make a series of center frequenies equally spaced in ERB-rate.
% 
% SYNTAX
% 
%     CFS = MAKEERBCFS(LOWCF,HIGHCF,NUMCHANS)
% 
% 
% DESCRIPTION
%
% CFS = MAKEERBCFS(LOWCF,HIGHCF,NUMCHANS) makes a vector of NUMCHANS center
% frequenies equally spaced on the ERB-rate scale between LOWCF and HIGHCF.
%
% Adapted from code written by: Guy Brown, University of Sheffield and 
% Martin Cooke

cfs = ErbRateToHz(linspace(HzToErbRate(mincf),HzToErbRate(maxcf),numchans));

% Copyright 2010, Chris Hummersone and University of Surrey
% 
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, 
% EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
% MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT 
% OF THIRD PARTY RIGHTS. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS 
% BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION 
% OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
% WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.