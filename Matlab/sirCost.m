function NLL = sirCost(tspan, params, data1, x0fcn, y2fcn, model)

    % Ensure parameters are non-negative
    params = abs(params);

    % ---------------------------------------------------------------
    % Solve model
    % ---------------------------------------------------------------
    [~, x] = ode45( model, tspan, [x0fcn(params), params], [], params);        
                              
    % ---------------------------------------------------------------
    % Model observations
    % ---------------------------------------------------------------
    y1 = y2fcn(x, params);

    % Avoid log(0) and numerical negative values
    y1 = max(y1, eps);

    % ---------------------------------------------------------------
    % Dimension check
    % ---------------------------------------------------------------
    if ~isequal(size(y1), size(data1))
        error('sirCost:DimensionMismatch', ...
            'Model output and observed data must have the same dimensions.');
    end

    % ---------------------------------------------------------------
    % Poisson negative log-likelihood
    %
    % NLL = sum(lambda - y*log(lambda))
    %
    % log(y!) is omitted because it does not depend on parameters.
    % ---------------------------------------------------------------
    NLL = sum(y1 - data1 .* log(y1), 'all');

    % ---------------------------------------------------------------
    % Numerical safety
    % ---------------------------------------------------------------
    if ~isfinite(NLL)
        NLL = realmax;
    end

end