function v= Variance_Explained_SVD(model,data,sigma,opt,n_para,n_data)

% function v = Variance_Explained_SVD(model,data,sigma,opt,n_para,n_data,n_rep)
% real sigma (std error of mean) needed, not just 1/sqrt(n_data) for opt='explainable'!!
% opt: 'orthodox' (default)
% opt: 'explainable-' : estimator with unbiased numerator and denominator
%                       (not much better MSE than orthodox)
% opt: 'explainable+' : estimator with conditioned/regularized denominator
%                       (much better MSE but based purely on simulations so far)

if nargin<3, sigma=1; end
if nargin<4, opt='orthodox'; end

model=model(:);
data =data (:);
sigma=sigma(:);

if nargin<6 || isempty(n_data), n_data=numel(data); end

switch opt
  
  case 'orthodox'
    v=1-var(data-model)/var(data);

  case 'explainable-'
    if nargin<6, n_data=length(model); end
    if numel(sigma)>1
      warning(['sigma must be scalar! using mean of range: ' num2str([min(sigma(:)) max(sigma(:))])]);
      sigma=mean(sigma(:));
    end
    v=1-(var(data-model)-(n_data-n_para)/(n_data-1)*sigma^2)/...
      (var(data)-sigma^2);
    
  case 'explainable+'
    if nargin<6, n_data=length(model); end
    if numel(sigma)>1
      warning(['sigma must be scalar! using mean of range: ' num2str([min(sigma(:)) max(sigma(:))])]);
      sigma=mean(sigma(:));
    end
    v=1-(var(data-model)-(n_data-n_para)/(n_data-1)*sigma^2)/...
      (var(data)+n_data/(n_data-1)*sigma^2);

  otherwise
    error(['invalid option ' opt]);
end
