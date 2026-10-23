function [x_hat,y_hat,theta_hat,P,S,x_hat_i,ensemble_history,param_corr,diagnostics] = ...
    enkfnew2026(z,y2fcn,Q,V_pram,x0fcn,paramests,n,model,error1)
%ENKFNEW2026 Physically constrained perturbed-observation EnKF.
% Augmented coordinates: [N,S,I,Q,R,beta,beta_f,mu_f,I0].
% Observation 1 corresponds to the initial date (elapsed day zero).
    paramests = abs(paramests(:)');
    n_mes = size(z,2);
    n_obs = size(z,1);
    n_param = numel(paramests);
    n_aug = 5+n_param;
    if isempty(z) || any(~isfinite(z(:))) || error1 < 0
        error('enkfnew2026:Observations','Supply finite observations and nonnegative observation error.');
    end
    [x_hat_i,initial] = initialize_epidemic_ensemble(x0fcn,paramests,V_pram,n);
    diagnostics = epidemic_constraint_diagnostics(initial,n_mes);
    x_hat = zeros(n_aug,n_mes);
    y_hat = zeros(n_obs,n_mes);
    theta_hat = zeros(n_param,n_mes);
    P = zeros(n_aug,n_aug,n_mes);
    S = zeros(n_obs,n_obs,n_mes);
    ensemble_history = zeros(n_aug,n,n_mes);
    param_corr = zeros(n_param,n_param,n_mes);

    for t = 1:n_mes
        % No model advance before assimilating the initial-date observation.
        if t > 1
            if iscell(Q), process_cov = Q{t}; else, process_cov = Q; end
            w = mvnrnd(zeros(1,n_aug),process_cov,n)';
            for i = 1:n
                x_hat_i(:,i) = propagate_epidemic_state(model,[t-2,t-1], ...
                    x_hat_i(:,i),x_hat_i(6:end,i));
            end
            [x_hat_i,correction] = project_epidemic_state(x_hat_i+w);
            diagnostics.prediction_corrected_members(t) = correction.corrected_count;
        end

        yk = z(:,t);
        R_k = diag(max(error1*abs(yk),1).^2);
        v = mvnrnd(zeros(1,n_obs),R_k,n)';
        y_members = observe_members(x_hat_i,y2fcn,n_obs);
        x_bar = mean(x_hat_i,2);
        y_bar = mean(y_members,2);
        Ex = x_hat_i-x_bar;
        Ey = y_members-y_bar;
        P_yt = Ey*Ey'/(n-1)+R_k;
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
        [x_hat_i,correction] = project_epidemic_state( ...
            x_hat_i+gain*(yk+v-y_members));
        diagnostics.analysis_corrected_members(t) = correction.corrected_count;

        % Recompute all reported moments after physical projection.
        y_members = observe_members(x_hat_i,y2fcn,n_obs);
        x_hat(:,t) = project_epidemic_state(mean(x_hat_i,2));
        theta_hat(:,t) = mean(x_hat_i(6:end,:),2);
        y_hat(:,t) = mean(y_members,2);
        Ex = x_hat_i-mean(x_hat_i,2);
        Ey = y_members-y_hat(:,t);
        P(:,:,t) = Ex*Ex'/(n-1);
        S(:,:,t) = Ey*Ey'/(n-1);
        ensemble_history(:,:,t) = x_hat_i;
        Ctheta = corrcoef(x_hat_i(6:end,:)');
        Ctheta(~isfinite(Ctheta)) = 0;
        Ctheta(1:n_param+1:end) = 1;
        param_corr(:,:,t) = Ctheta;
        diagnostics.minimum_Sf(t) = min(x_hat_i(1,:)-sum(x_hat_i(2:5,:),1));
        diagnostics.minimum_population(t) = min(x_hat_i(1,:));
    end
end

function observations = observe_members(members,y2fcn,n_obs)
    observations = zeros(n_obs,size(members,2));
    for i = 1:size(members,2)
        observations(:,i) = y2fcn(members(1:5,i)',members(6:end,i)');
    end
end