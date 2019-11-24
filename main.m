function TestAdCampaign
% Update info:
% 1. changed: n, abar,cs, ds, ps, B, ahat
%
%This file requires YALMIP to be installed
%https://yalmip.github.io/

n = 3;
abar = [0.86 0.39 0.38]';
cs = 18525*ones(n,1);
ds = 15257*ones(n,1);
ps = ones(n,1);
B = 1e7;
rho = -1;
ahat = ([0.54 0.11 0.08].*abar')';
Gamma = 1;

x0 = B/n*ones(n,1); 

%solve nominal problem
[fval,x] = solveNominalProb(abar,cs,ds,ps,B,[]);
nomSol.x = x; nomSol.fval = fval;

%solve robust formulation under Box uncertainty set
[fval,x] = solveNominalProb(abar-ahat,cs,ds,ps,B,[]);
boxSol.x = x; boxSol.fval = fval;

%robust robust formulation under Budgeted uncertainty set
[fval,x] = solveRobustProbBudgetSet(abar,ahat,Gamma,cs,ds,ps,B,[]);
budgetSol.x = x; budgetSol.fval = fval;

%evaluate nominal values
nomSol.nomVal = nomSol.fval;
boxSol.nomVal = solveNominalProb(abar,cs,ds,ps,B,boxSol.x);
budgetSol.nomVal = solveNominalProb(abar,cs,ds,ps,B,budgetSol.x);

%evaluate robust values under budgeted uncertainty set
nomSol.wcVal = solveRobustProbBudgetSet(abar,ahat,Gamma,cs,ds,ps,B,nomSol.x)
boxSol.wcVal = solveRobustProbBudgetSet(abar,ahat,Gamma,cs,ds,ps,B,boxSol.x)
budgetSol.wcVal = budgetSol.fval;

display('These are the three solutions founds : "nominal"  "robust with box" "robust with budget"');
[nomSol.x boxSol.x budgetSol.x]

display(sprintf('The nominal solution predicts a profit of %0.3f, its nominal profit is %0.3f, its worst-case profit is %0.3f',nomSol.fval,nomSol.nomVal,nomSol.wcVal));
display(sprintf('The robust with box solution predicts a profit of %0.3f, its nominal profit is %0.3f, its worst-case profit is %0.3f',boxSol.fval,boxSol.nomVal,boxSol.wcVal));
display(sprintf('The robust with budget solution predicts a profit of %0.3f, its nominal profit is %0.3f, its worst-case profit is %0.3f',budgetSol.fval,budgetSol.nomVal,budgetSol.wcVal));

return

function [fval,x] = solveNominalProb(abar,cs,ds,ps,B,xfixed)
n = length(abar);
%nominal problem
t = sdpvar(1,1);
y = sdpvar(n,1);
x = sdpvar(n,1);
tmp = sdpvar(n,1);
objective = -t;
constraints = [y == 1+x./ds, ps'*x<= B, x>=0];
constraints = constraints + [x >= 0.1];
for k=1:n
    constraints = constraints + [tmp(k) <= cs(k)*cpower(y(k),abar(k))-cs(k)];
end;
constraints = constraints + [t <= sum(tmp)];
if ~isempty(xfixed), constraints = constraints + [x == xfixed]; end;
options = sdpsettings('allownonconvex',0);
sol = optimize(constraints, objective, options);
x = value(x);
fval = value(t);
return


function [fval,x] = solveRobustProbBudgetSet(abar,ahat,Gamma,cs,ds,ps,B,xfixed)
n = length(abar);
t = sdpvar(1,1);
x = sdpvar(n,1);
s = sdpvar(1,1);
v = sdpvar(n,1);
w = sdpvar(n,1);
lambda = sdpvar(1,1);
tmp2 = sdpvar(n,1);
objective = -t;
constraints = [ps'*x<= B, x>=0];
constraints = constraints + [x >= 0.1];
constraints = constraints + [v<=0, lambda>= -diag(ahat)*v-w, lambda >=0, w>=0];
constraints = constraints + [s>=Gamma*lambda+sum(w)];
for k=1:n
constraints = constraints + [tmp2(k)>=cs(k)*fancyFunction([-v(k)/cs(k); x(k)./ds(k)])+cs(k)];
end;
constraints = constraints + [t+abar'*v+s+sum(tmp2)<=0];
if ~isempty(xfixed), constraints = constraints + [x == xfixed]; end;
options = sdpsettings('allownonconvex',0);
options.fmincon.MaxIter = 4000;
options.fmincon.TolCon = 1.0000e-12;
options.fmincon.TolFun = 1.0000e-12;
options.fmincon.TolX = 1.0000e-16;
options.fmincon.MaxFunEvals = 8000;
options.fmincon.Algorithm = 'trust-region-reflective';
sol = optimize(constraints, objective, options);
x = value(x);
fval = value(t);
return


