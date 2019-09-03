function [x_keep, logL_keep, logQ_keep, accrate] = mcmc(func,D,X,Niter,verbose)

% TODO:
%  * Help description
%  * Option for computing prediction at each iteration?

% Parse out data and errors
d       = D.d;
Sig     = D.Sig;

% Parse out model initiual guess, step size, bounds and priors
x0      = X.x0;
xstep   = X.xstep;
xbnds   = X.xbnds;
xprior  = X.xprior;
C       = X.C;

% Make sure initial guess is within bounds
if ~isempty(xbnds)
    assert(all(x0 > xbnds(:,1) & x0 < xbnds(:,2)),'Initial guess is outside bounds')
end

% Initialize
x       = x0;
count   = 0;
x_keep  = cell(Niter,1);
logL_keep  = nan(Niter,1);
logQ_keep  = nan(Niter,1);

fun = fcnchk(func);

% Check if likelihood is uniform
if isempty(Sig) | isinf(Sig)
    logLc = 0;
    logLx = 0;
    compute_likelihood = 0;
    if verbose, disp('Likelihood is uniform'), end
else
    compute_likelihood = 1;
    if verbose, disp('Likelihood is not uniform'), end
end

% Check if prior is non-informative
if isempty(xprior) || isempty(C) || all(all(isinf(C)))
    logPc = 0;
    logPx = 0;
    compute_prior = 0;  % no need to compute prior since logP = 0
    if verbose, disp('Prior is uniform for all parameters'), end
else
    compute_prior = 1;
    if verbose, disp('Prior is not uniform for all parameters'), end
end

[f1,f2,f3,f4,f5,f6,f7,f8]=deal(1);

if verbose, disp('Sampling started...'), end
for i = 1:Niter

    % Draw candidate model from proposal distribution
    if 0
        c = x + xstep.*randn(length(x),1);      % Gaussian distribution
    else
        c = x + xstep.*(rand(length(x),1)-0.5); % Uniform distribution
    end

    % Check bounds
    if ~isempty(xbnds) & and(c > xbnds(:,1), c < xbnds(:,2))

        % Compute predicted data
        Gc = fun(c);   % Predicted data of proposed model
        Gx = fun(x);   % Predicted data of current model

        % Compute log likelihood
        if compute_likelihood
            if numel(Sig) == 1
                if verbose & f1==1, disp('Single (scalar) error provided.'), f1 = 0; end
                % scalar data error
                logLc = -0.5*(norm(d-Gc))^2/Sig^2;  % proposed model
                logLx = -0.5*(norm(d-Gx))^2/Sig^2;  % current model
            else
                if verbose & f2==1, disp('Data covariance matrix provided.'), f2 = 0; end
                % data covariance matrix
                rc = d - Gc;               % residual from proposed model
                logLc = -0.5*rc'*(Sig\rc);  % proposed model

                rx = d - Gx;               % residual from current model
                logLx = -0.5*rx'*(Sig\rx);  % proposed model
            end
        end

        % Compute log prior
        if compute_prior

            k = numel(xprior);

            % Mixed prior (first k parameters have Gaussian prior, rest
            % have uniform)
            if k < numel(x0)
                if verbose & f3==1
                    disp(['Mixed prior, with first ' num2str(k) ' parameters having Gaussian priors.'])
                    f3 = 0;
                end

                if numel(C) == 1
                    % scalar model error
                    if verbose & f4==1, disp('Single (scalar) model error provided.'), f4 = 0; end
                    logPc = -0.5*(norm(xprior(1:k) - c(1:k)))^2/C^2;  % proposed model
                    logPx = -0.5*(norm(xprior(1:k) - x(1:k)))^2/C^2;  % current model

                else
                    % model covariance matrix
                    if verbose & f5==1, disp('Model covariance matrix provided.'), f5 = 0; end
                    mc = c(1:k) - xprior(1:k); % residual from proposed model
                    logPc = -0.5*mc'*(C\mc);  % proposed model

                    mx = x(1:k) - xprior(1:k);  % residual from current model
                    logPx = -0.5*mx'*(C\mx);  % proposed model
                end


            else
                % Gaussian prior for all parameters
                if verbose & f6==1
                    disp('Gaussian prior for all parameters')
                    f6 = 0;
                end

                if numel(C) == 1
                    % scalar model error
                    if verbose & f7==1, disp('Single (scalar) model error provided.'), f7=0; end
                    logPc = -0.5*(norm(xprior-c))^2/C^2;  % proposed model
                    logPx = -0.5*(norm(xprior-x))^2/C^2;  % current model
                else
                    % model covariance matrix
                    if verbose & f8==1, disp('Model covariance matrix provided.'), f8=0; end
                    mc = c - xprior;               % residual from proposed model
                    logPc = -0.5*mc'*(C\mc);  % proposed model

                    mx = x - xprior;               % residual from current model
                    logPx = -0.5*mx'*(C\mx);  % proposed model
                end
            end
        end


        % Compute log posteriors
        logQc = logLc + logPc; % proposed model
        logQx = logLx + logPx; % current model

        % Acceptance step
        postratio = exp(logQc - logQx);  % posterior ratio
        Pacc = min(postratio, 1);

        % Acceptance step
        if Pacc >= rand
            x_keep{i} = c;          % Keep proposed model
            logL_keep(i) = logLc;   % Log likelihood of accepted model
            logQ_keep(i) = logQc;   % Log posterior of accepted model
            count = count+1;        % Update count of accepted models
        else
            x_keep{i} = x;          % Reject proposed model, keep current model
            logL_keep(i) = logLx;   % Log likelihood of current model
            logQ_keep(i) = logQx;   % Log posterior for current model
        end

    else

        x_keep{i} = x;          % Reject proposed model, keep current model
        logL_keep(i) = logLx;
        logQ_keep(i) = logQx;

    end

    x = x_keep{i};          % Update current model

    pct = i/Niter*100;
    if verbose & mod(pct,10)==0
        disp([num2str(pct) '% done. Cumulative acceptance rate = ' ...
            num2str(100*count/i) '%'])
    end

end

if verbose, disp('Sampling ended'), end

accrate = count/Niter;  % acceptance ratio

if verbose, disp(['Acceptance rate = ',num2str(accrate*100),' %']), end
