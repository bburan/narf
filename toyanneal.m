% Search vector tests 


%  2. Samples according to a gaussian with mean mu and covariance matrix Sig
sample_gaussianly = @(mu, Sig, n_samps) (repmat(mu,n_samps,1) + randn(n_samps,2)*chol(Sig));

% Single Constraint At A Time 
% sample_scaat = @();

% Generate values from a bivariate normal distribution with specified mean
% vector and covariance matrix.
mu = [10 20];
Sigma = [1 .8; 
         .8 2];

samps_u = sample_uniformly([-5 1], [5 3], 1000);
samps_g = sample_gaussianly(mu, Sigma, 1000);
     
figure;
plot(samps_u(:,1), samps_u(:,2), 'b.', ...
    samps_g(:,1), samps_g(:,2), 'g.');

% --------------------------------
% Full fledged example

% Hidden information that we don't know
hidden = [1.2 2.4];
evalfn = @(guess) sqrt(sum((guess - hidden).^2));

% A demonstration of how to find the hidden information
n_samps = 1000;
best = [0 0];
bestscore = inf;

% Uniform sample algorithm
testpoints = sample_uniformly([-5 1], [5 3], n_samps);
for i = 1:n_samps;
    p = testpoints(i,:);
    score = evalfn(p);
    if (score < bestscore)
        bestscore = score;
        best = p;
    end
end

disp('Best value found uniformly:'); disp(best);

% Second stage, gaussian sampling near best point
testpoints = sample_gaussianly(best, [.1 0.01; 0.01 .1], n_samps);
for i = 1:n_samps;
    p = testpoints(i,:);
    score = evalfn(p);
    if (score < bestscore)
        bestscore = score;
        best = p;
    end
end

disp('After Gaussian 2nd stage:'); disp(best);

% Third stage, a stepping algorithm
[x_bst, s_bst] = scaat_stepper([0 0], evalfn([0 0]), evalfn, 20, 1);

disp('Scaat stepper:'); disp(x_bst);


