function [strfH, strfHJN, strfHJN_std] = cal_Strf(fstim, fstim_spike,...
     stim_spike_JNf, stim_size, stim_spike_size, stim_spike_JNsize,...
     nb, nt, nJN, tol, save_flag)
%
%  [strfH, strfHJN, strfHJN_std] = cal_Strf(fstim, fstim_spike,...
%     stim_spike_JNf, stim_size, stim_spike_size, stim_spike_JNsize,...
%     nb, nt, nJN, tol)
%      -- Calculate strf, JackKnifed strf and contour JN strf
%     Input:
%         fstim: FFT of stim auto-correlation
%         fstim_spike: FFT of stim_spike cross-correlation
%         stim_spike_JNf: FFT of JackKnifed stim_spike cross-correlation
%         stim_size: size of stim auto-correlation
%         stim_spike_size: size of stim_spike cross-correlation
%         stim_spike_JNsize: size of JN stim_spike cross-correlation 
%         nb: length of spatio domain
%         nt: length of time domain
%         nJN: num of JackKnife case
%         tol: one tol. value
%         save_flag: the flag to save the intermediate result for specific tol
%                    The default value is 0 (no save), 1 otherwise.
%     Output:
%         strfH: the estimated strf
%         strfHJN: the estimated JackKnifed version strf
%         strfHJN_std: the estimated JN strf
%         
%             STRFPAK: STRF Estimation Software
% Copyright ©2003. The Regents of the University of California (Regents).
% All Rights Reserved.
% Created by Theunissen Lab and Gallant Lab, Department of Psychology, Un
% -iversity of California, Berkeley.
%
% Permission to use, copy, and modify this software and its documentation
% for educational, research, and not-for-profit purposes, without fee and
% without a signed licensing agreement, is hereby granted, provided that
% the above copyright notice, this paragraph and the following two paragr
% -aphs appear in all copies and modifications. Contact The Office of Tec
% -hnology Licensing, UC Berkeley, 2150 Shattuck Avenue, Suite 510,
% Berkeley, CA 94720-1620, (510) 643-7201, for commercial licensing
% opportunities.
%
%IN NO EVENT SHALL REGENTS BE LIABLE TO ANY PARTY FOR DIRECT, INDIRECT,
%SPECIAL, INCIDENTAL, OR CONSEQUENTIAL DAMAGES, INCLUDING LOST PROFITS,
%ARISING OUT OF THE USE OF THIS SOFTWARE AND ITS DOCUMENTATION, EVEN IF
%REGENTS HAS BEEN ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
%
%REGENTS SPECIFICALLY DISCLAIMS ANY WARRANTIES, INCLUDING, BUT NOT
%LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
%PARTICULAR PURPOSE. THE SOFTWARE AND ACCOMPANYING DOCUMENTATION, IF ANY,
%PROVIDED HEREUNDER IS PROVIDED "AS IS". REGENTS HAS NO OBLIGATION TO PRO
%-VIDE MAINTENANCE, SUPPORT, UPDATES, ENHANCEMENTS, OR MODIFICATIONS.



% =======================================================
% Check if we have all input
% =======================================================
if ~exist('save_flag')
    save_flag = 0;
end

% ========================================================
% Forward Filter - The algorithm is from FET's filters2.m
% ========================================================
nf = (nt-1)/2 +1;

% Allocate space for all arrays
clear j;
stim_mat = zeros(nb,nb)+j.*zeros(nb,nb);
cross_vect = zeros(nb,1);
cross_vectJN = zeros(nJN,nb);
h = zeros(1,nb);
hJN = zeros(nJN, nb);
ffor=zeros(stim_spike_size);
fforJN=zeros(stim_spike_JNsize);
strfH =zeros(stim_spike_size);
strfHJN=zeros(stim_spike_JNsize);
cums=zeros(nf,nb+1);
ranktest=zeros(1,nf);
stimnorm=zeros(1,nf);
clear j;


% ========================================================
% Find the maximum norm of all the matrices
% ========================================================
for iff=1:nf
    nc = 1;
    for ii=1:nb
	for jj=ii:nb
        stim_mat(ii,jj) = fstim(nc,iff);
	    if ii ~= jj
		stim_mat(jj,ii) = conj(fstim(nc,iff));
	    end
	    nc = nc +1;
	end
    end
    stimnorm(1,iff)=norm(stim_mat);
end

% One way of doing it
ranktol=tol*max(stimnorm);

% Do the matrix inversion for each frequency 
for iff=1:nf
    % Stuff stim matrix and cross-correlation vectors
    nc = 1;
    for ii=1:nb
	for jj=ii:nb
            stim_mat(ii,jj) = fstim(nc,iff);
	    if ii ~= jj
		stim_mat(jj,ii) = conj(fstim(nc,iff));
	    end
	    nc = nc +1;
	end
	% cross_vect(ii) = conj(fstim_spike(ii,iff));
	cross_vect(ii) = fstim_spike(ii,iff);
	for iJN=1:nJN
            % cross_vectJN(iJN,ii) = conj(stim_spike_JNf(ii,iff,iJN));
	    cross_vectJN(iJN,ii) = stim_spike_JNf(ii,iff,iJN);
	end
    end

    % do an svd decomposition
    ranktest(1,iff)=rank(stim_mat,ranktol);
    [u,s,v] = svd(stim_mat);
    tots = s(1,1);
    cums(iff,2) = s(1,1);
    for ii=2:nb
        tots = tots + s(ii,ii);
        cums(iff,ii+1) = cums(iff,ii) + s(ii,ii);
    end
    is = zeros(nb,nb);
    for ii=1:nb+1
	cums(iff,ii) = cums(iff,ii)/tots;
    end

    % ncutt_off = round(1.0/cums(iff,2))+1;
    for ii=1:nb
        if ii>ranktest(1,iff)
            is(ii,ii)=(1.0/ranktol)*exp(-(ii-ranktest(1,iff))^2/8);
	else
	    is(ii,ii)=1.0/s(ii,ii);
	end
    end

    h = v*is*(u'*cross_vect);
    for  iJN=1:nJN
	hJN(iJN,:) = (v*is*(u'*cross_vectJN(iJN,:).')).';
    end
   
    for ii=1:nb
	ffor(ii,iff) = h(ii);
	fforJN(ii,iff,:) = hJN(:,ii);
	if iff ~= 1
            ffor(ii,nt+2-iff) = conj(h(ii));
	    fforJN(ii,nt+2-iff,:) = conj(hJN(:,ii));
	end
    end
end

for ii=1:nb
   strfH(ii,:) = real(ifft(ffor(ii,:)));
   for iJN=1:nJN
       strfHJN(ii,:,iJN) = real(ifft(fforJN(ii,:,iJN)));
   end
end

strfHJN_mean = mean(strfHJN,3);
strfHJN_var = zeros(size(strfHJN_mean));
strfHJN_nJN = zeros(size(strfHJN_mean));

for iJN=1:nJN
    strfHJN_var = strfHJN_var + (strfHJN(:,:,iJN) - strfHJN_mean).^2;
end

strfHJN_var = (1-1/nJN)*strfHJN_var;
strfHJN_std = strfHJN_var.^.5;

% =======================================================
% Save the result into Output/files
% =======================================================
if save_flag == 1
    currentPath = pwd;

    global outputPath
     if ~isempty(outputPath)
         cd (outputPath);
     else
         disp('Saving output to Output Dir.');
         stat = mkdir('Output');
         cd('Output');
         outputPath = pwd;
     end

    save('strfH.mat', 'strfH');
    save('strfH_std.mat', 'strfHJN_std');
    for iJN=1:nJN
	filename = sprintf('strfHJN%d.mat',iJN);
	strfHJN_nJN = strfHJN(:,:,iJN);
	save(filename, 'strfHJN_nJN');
    end
    cd(currentPath);
end

% =======================================================
% END OF CAL_STRF
% =======================================================
