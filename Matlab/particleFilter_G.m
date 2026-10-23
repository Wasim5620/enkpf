function [x_hatP, y_hatP, theta_hatP, ESS_history, unique_history, diagnostics] = particleFilter_G( ...
    z, y2fcn, Q, V_pram, x0fcn, paramests, n, model, error1)
%% ==========================================================
%  PF-G : Gaussian Particle Filter (Simon bootstrap / SIR)
% ==========================================================
    paramests   = abs(paramests);
    x0          = x0fcn(paramests);
    n_mes       = size(z,2);
    n_sta       = length(x0);
    n_obs_var   = size(z,1);
    n_obs_param = numel(paramests);
    n_sta_var   = n_sta + n_obs_param;
    
    x_hatP         = zeros(n_sta_var, n_mes);
    y_hatP         = zeros(n_obs_var, n_mes);
    theta_hatP     = zeros(n_obs_param, n_mes);
    ESS_history    = zeros(1, n_mes);
    unique_history = zeros(1, n_mes);
    
    % Couple each sampled I0 to its physical initial state.
    [x_hat_i, initial] = initialize_epidemic_ensemble(x0fcn, paramests, V_pram, n);
    x_hat_i(1:n_sta, :) = max(x_hat_i(1:n_sta, :) + 0.01*abs(x0(:)).*randn(n_sta,n), 0);
    [x_hat_i, initial_noise] = project_epidemic_state(x_hat_i);
    x_hat_i(9,:) = x_hat_i(3,:); % I0 must describe the noisy feasible initial I.
    diagnostics = epidemic_constraint_diagnostics(initial, n_mes);
    diagnostics.initial_corrected_members = ...
        diagnostics.initial_corrected_members + initial_noise.corrected_count;

    wts = ones(1,n)/n;
    for t = 1:n_mes
        if iscell(Q)
            processCov = Q{t};
        else
            processCov = Q;
        end
        if t > 1
            w = mvnrnd(zeros(n_sta_var,1), processCov, n)';
        else
            w = zeros(n_sta_var,n); % Initial-date observation: no advance.
        end
        yk  = z(:,t);
        R_k = construct_R_matrix(yk, n_obs_var, error1);
        R_k = (R_k + R_k')/2;
        
        % 2a. Time update / propagation
        y_hat_i = zeros(n_obs_var, n);
        for i = 1:n
            if t > 1
                x_hat_i(:,i) = propagate_epidemic_state(model, [t-2 t-1], ...
                    x_hat_i(:,i), x_hat_i(n_sta+1:end,i));
            end
            [x_hat_i(:,i), correction] = project_epidemic_state(x_hat_i(:,i) + w(:,i));
            diagnostics.prediction_corrected_members(t) = ...
                diagnostics.prediction_corrected_members(t) + correction.corrected_count;
            y_hat_i(:,i) = y2fcn(x_hat_i(1:n_sta,i)', x_hat_i(n_sta+1:end,i)');
        end
        
        % 2b. Gaussian likelihood
        loglik = zeros(1,n);
        [L,p] = chol(R_k,'lower');
        if p == 0
            logdetR = 2*sum(log(diag(L)));
            for i = 1:n
                v = L \ (yk - y_hat_i(:,i));
                loglik(i) = -0.5*(v'*v + logdetR + n_obs_var*log(2*pi));
            end
        else
            Rinv = pinv(R_k);
            logdetR = log(max(det(R_k),eps));
            for i = 1:n
                innov = yk - y_hat_i(:,i);
                loglik(i) = -0.5*(innov'*Rinv*innov + logdetR + n_obs_var*log(2*pi));
            end
        end
        
        % Normalize weights
        if any(~isfinite(loglik))
            loglik(~isfinite(loglik)) = -realmax;
        end
        % Retain the previous normalized weights when no resampling occurred.
        logw = log(max(wts,realmin)) + loglik;
        logw = logw - max(logw);
        wts = exp(logw);
        if sum(wts) <= 0 || any(~isfinite(wts))
            wts = ones(1,n)/n;
        else
            wts = wts/sum(wts);
        end
        
        % 2c. Statistics before resampling
        x_hatP(:,t)     = project_epidemic_state(x_hat_i * wts');
        y_hatP(:,t)     = y_hat_i * wts';
        theta_hatP(:,t) = x_hat_i(n_sta+1:end,:) * wts';
        
        ESS = 1/sum(wts.^2);
        ESS_history(t) = ESS;
        
        % 2d. Resample + mild roughening
        if ESS < 0.5*n
            [x_hat_i, n_uniq] = systematicResample(x_hat_i, wts);
            unique_history(t) = n_uniq;
            
            Krough = 0.3;
            for d = 1:n_sta_var
                range_d = max(x_hat_i(d,:)) - min(x_hat_i(d,:));
                sigma_d = Krough * range_d * n^(-1/n_sta_var);
                x_hat_i(d,:) = x_hat_i(d,:) + sigma_d*randn(1,n);
            end
            [x_hat_i, correction] = project_epidemic_state(x_hat_i);
            diagnostics.resampling_corrected_members(t) = correction.corrected_count;
            wts = ones(1,n)/n; % Resampled particles start with equal weights.
        else
            unique_history(t) = n; % No resampling needed, all particles unique
        end
        diagnostics.minimum_Sf(t) = min(x_hat_i(1,:) - sum(x_hat_i(2:5,:),1));
        diagnostics.minimum_population(t) = min(x_hat_i(1,:));
    end
end

%% Systematic resampling with unique count return
function [xr, num_unique] = systematicResample(X, w)
    N = size(X,2);
    edges = [0, cumsum(w)]; edges(end) = 1;
    u0 = rand/N; idx = zeros(1,N); j = 1;
    for i = 1:N
        ui = u0 + (i-1)/N;
        while ui > edges(j+1), j = j+1; end
        idx(i) = j;
    end
    num_unique = length(unique(idx));
    xr = X(:,idx);
end

%% Relative Gaussian observation covariance
function R = construct_R_matrix(yk, nobs, err)
    if nobs == 1
        R = max(err*abs(yk(1)),1)^2;
    elseif nobs == 2
        R = diag([max(err*abs(yk(1)),1)^2, max(err*abs(yk(2)),1)^2]);
    elseif nobs == 3
        R = diag([max(err*abs(yk(1)),1)^2, ...
                  max(err*abs(yk(2)),1)^2, ...
                  max(err*abs(yk(3)),1)^2]);
    else
        error('Only 1-3 observation variables supported.');
    end
end
