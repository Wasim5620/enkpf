%cost function
function P_theta = sirCost_LSM(tspan,params,data1,x0fcn,y2fcn, model)
params = abs(params);
[~,x] = ode45(model,tspan,[x0fcn(params), params],[],params);
y1 = y2fcn(x,params);

merr = (y1-data1);
P_theta = sum(merr.^2, "all");

