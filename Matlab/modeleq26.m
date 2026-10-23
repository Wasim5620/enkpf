function dxdt = modeleq26(~, x, params)
%% Known (fixed) parameters
Lambda = (5.9/1000)*(1/365);
mu     = (5.7/1000)*(1/365);
gamma  = 1/4;
alpha  = 1/14;
mmu    = 0.0221;          % disease-induced mortality
delta  = 1/1000;

%% Parameters to estimate
beta  = params(1);
betaf = params(2);
muf   = params(3);
% rbeta = 0.02;
rbeta = params(4);        % <-- now estimated

% Enforce the constraint 0 ≤ rbeta ≤ 1
rbeta = max(0.001, min(1, rbeta));

%% State variables
N  = x(1);
S  = x(2);
I  = x(3);
Q  = x(4);
R  = x(5);
Sf = N - (S + I + Q + R);

%% ODEs
dxdt = zeros(numel(x),1);

dxdt(1) = Lambda*N - mu*N - mmu*Q;

dxdt(2) = Lambda*N ...
    - (beta*S*I / N) ...
    - betaf*S*(1 - exp(-delta*Q)) ...
    + (muf*Sf*(S+R)/N) ...
    - mu*S;

dxdt(3) = (beta*S*I / N) + (rbeta*beta*Sf*I / N) - (gamma + mu)*I;

dxdt(4) = gamma*I - (mu + alpha + mmu)*Q;

dxdt(5) = alpha*Q - mu*R;

% Dummy equations for the parameters (if you are using an augmented state)
    % if numel(x) > 5
    %     dxdt = [dxdt; zeros(numel(y)-5,1)];
    % end
end