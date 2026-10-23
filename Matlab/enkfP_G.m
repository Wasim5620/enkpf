function [x_hat,y_hat,theta_hat,P,S,x_hat_i_final,theta_hat_i_final,ens_his,param_corr,diagnostics] = ...
    enkfP_G(z,y2fcn,Q,V_pram,x0fcn,paramests,n,model,error1)
%ENKFP_G Constrained project-specific EnKF/residual-weighted hybrid filter.
% Retains the existing residual-score/resampling scheme. This is not the
% conventional tempered-likelihood Frei-Kuensch EnKPF implementation.
% Augmented coordinates: [N,S,I,Q,R,beta,beta_f,mu_f,I0].
    paramests = abs(paramests(:)');
    n_mes = size(z,2);
    n_obs = size(z,1);
    n_param = numel(paramests);
    n_aug = 5+n_param;
    if isempty(z) || any(~isfinite(z(:))) || error1 < 0
        error('enkfP_G:Observations','Supply finite observations and nonnegative observation error.');
    end
    [members,initial] = initialize_epidemic_ensemble(x0fcn,paramests,V_pram,n);
    diagnostics = epidemic_constraint_diagnostics(initial,n_mes);
    diagnostics.ESS = zeros(1,n_mes);
    diagnostics.unique_parents = n*ones(1,n_mes);
    diagnostics.weight_underflow_uniform_fallback = false(1,n_mes);
    x_hat = zeros(n_aug,n_mes);
    y_hat = zeros(n_obs,n_mes);
    theta_hat = zeros(n_param,n_mes);
    P = zeros(n_aug,n_aug,n_mes);
    S = zeros(n_obs,n_obs,n_mes);
    ens_his = zeros(n_aug,n,n_mes);
    param_corr = zeros(n_param,n_param,n_mes);

    for t = 1:n_mes
        if t > 1
            if iscell(Q), process_cov = Q{t}; else, process_cov = Q; end
            w = mvnrnd(zeros(1,n_aug),process_cov,n)';
            for i = 1:n
                members(:,i) = propagate_epidemic_state(model,[t-2,t-1], ...
                    members(:,i),members(6:end,i));
            end
            [members,correction] = project_epidemic_state(members+w);
            diagnostics.prediction_corrected_members(t) = correction.corrected_count;
        end
        % Preserve the existing parameter variance floor, then re-constrain.
        if trace(cov(members(6:end,:)')) < 1e-6
            members(6:end,:) = members(6:end,:)+1e-4*randn(n_param,n);
            [members,correction] = project_epidemic_state(members);
            diagnostics.prediction_corrected_members(t) = ...
                diagnostics.prediction_corrected_members(t)+correction.corrected_count;
        end

        yk = z(:,t);
        R_k = diag(max(error1*abs(yk),1).^2);
        v = mvnrnd(zeros(1,n_obs),R_k,n)';
        R_ek = v*v'/(n-1);
        y_members = observe_members(members,y2fcn,n_obs);
        x_bar = mean(members,2);
        y_bar = mean(y_members,2);
        Ex = members-x_bar;
        Ey = y_members-y_bar;
        P_yt = Ey*Ey'/(n-1)+R_ek;
        P_yt = (P_yt+P_yt')/2;
        if rcond(P_yt) < 1e-12
            P_yt = P_yt+max(1e-8,1e-6*norm(P_yt,'fro'))*eye(n_obs);
        end
        P_xyt = Ex*Ey'/(n-1);
        [Rchol,flag] = chol(P_yt);
        if flag == 0
            gain = (P_xyt/Rchol)/Rchol';
        else
            gain = P_xyt*pinv(P_yt,1e-12);
        end
        [members,correction] = project_epidemic_state( ...
            members+gain*(yk+v-y_members));
        diagnostics.analysis_corrected_members(t) = correction.corrected_count;
        y_members = observe_members(members,y2fcn,n_obs);
        residuals = yk-y_members;

        % Original Gaussian residual score, accompanied by diagnostics.
        scale = std(residuals(:))/n_aug;
        if ~isfinite(scale) || scale <= 0, scale = 1e-6; end
        weights = normpdf(vecnorm(residuals,2,1),0,scale);
        if ~isfinite(sum(weights)) || sum(weights) <= 0
            weights = ones(1,n)/n;
            diagnostics.weight_underflow_uniform_fallback(t) = true;
        else
            weights = weights/sum(weights);
        end
        diagnostics.ESS(t) = 1/sum(weights.^2);
        if diagnostics.ESS(t) < 0.5*n
            param_cov = cov(members(6:end,:)');
            param_cov = (param_cov+param_cov')/2;
            edges = [0,cumsum(weights)];
            edges(end) = 1.0000001;
            indices = discretize((rand()+(0:n-1))/n,edges);
            diagnostics.unique_parents(t) = numel(unique(indices));
            members = members(:,indices);
            jitter = mvnrnd(zeros(1,n_param), ...
                0.05*param_cov+1e-10*eye(n_param),n)';
            members(6:end,:) = members(6:end,:)+jitter;
        end
        [members,correction] = project_epidemic_state(members);
        diagnostics.resampling_corrected_members(t) = correction.corrected_count;

        % Moments describe the final constrained/resampled ensemble.
        y_members = observe_members(members,y2fcn,n_obs);
        x_hat(:,t) = project_epidemic_state(mean(members,2));
        theta_hat(:,t) = mean(members(6:end,:),2);
        y_hat(:,t) = mean(y_members,2);
        Ex = members-mean(members,2);
        Ey = y_members-y_hat(:,t);
        P(:,:,t) = Ex*Ex'/(n-1);
        S(:,:,t) = Ey*Ey'/(n-1);
        ens_his(:,:,t) = members;
        Ctheta = corrcoef(members(6:end,:)');
        Ctheta(~isfinite(Ctheta)) = 0;
        Ctheta(1:n_param+1:end) = 1;
        param_corr(:,:,t) = Ctheta;
        diagnostics.minimum_Sf(t) = min(members(1,:)-sum(members(2:5,:),1));
        diagnostics.minimum_population(t) = min(members(1,:));
    end
    x_hat_i_final = members;
    theta_hat_i_final = members(6:end,:);
end

function observations = observe_members(members,y2fcn,n_obs)
    observations = zeros(n_obs,size(members,2));
    for i = 1:size(members,2)
        observations(:,i) = y2fcn(members(1:5,i)',members(6:end,i)');
    end
end