function FIM = MiniFisher1(T, paramests, x0fcn, xfcn,y1fcn)

delta = 0.001; % Percent change we'll use, default 0.1%

% Simulate the system to determine the size of the state vector
 [~, x1] = ode45(xfcn, T, [x0fcn(paramests), paramests], [], paramests);

% Determine the number of states and the length of the time span
% num_states = size(x_test, 2);
num_timepoints = length(T);

% Preallocate the design matrix X
% Assuming we're interested in the sensitivities of the 7th and 8th states
num_sensitivities = size(y1fcn(x1,paramests),2);
X = zeros(num_timepoints * num_sensitivities, length(paramests));

for j = 1:length(paramests)
    params1 = paramests;
    params2 = paramests;
    params1(j) = (1 + delta) * paramests(j);
    params2(j) = (1 - delta) * paramests(j);
    
    [~, x1] = ode45(xfcn, T, [x0fcn(params1), params1], [], params1);
    [~, x2] = ode45(xfcn, T, [x0fcn(params2), params2], [], params2);
    
    % Calculate the sensitivities for the 7th and 8th states
    Y = (y1fcn(x1,params1) - y1fcn(x2,params2)) / (2 * delta * paramests(j));
    
    % Fill in the jth column of the design matrix with the sensitivities to parameter j
    X(:, j) = reshape(Y, [], 1);
end

% FIM (simplified w/o weighting term)
FIM = X' * X;

end
