function Judge = bowlf(Step, StepSize, Data0, Data, flag)
%BOSTWLF Judge whether Data satisfy the Wolfe Criterion and the basic
%criterion in inexact line search
%
% In line search problem P: a = argmin f(x+ad), along with g the gradient
% of f, the basic criterion is
%       F_new <= F + rho*a*g'*d,    rho in (0,.5)   (C1)
% and the Wolfe Criterion is
%       g_new'*d >= sigma*g'*d,     sigma in (0,1)  (C2)
% Notice that: sigma > rho
%
% Call
% Judge = bowlf(Step, StepSize, Data0, Data)
% Judge = bowlf(Step, StepSize, Data0, Data, flag)
%
% Input
% Step:     d in P, n-vector
% StepSize: a in P, scalar
% Data:     F & g in P & C1 & C2, struct with fields F g
% New:      F_new & g_new in P & C1 & C2, struct with fields F g
% flag:     rho & sigma in C1 & C2, array
%       flag(1):    sigma in C2
%           (2):    rho in C1
%       default:    [0.95 0.05]
% Output
% Judge:    answer to judge the satisfication, logical number 0/1
%       =1: just satisfy the criterions
%       =0: fail to satisfy those
%

% Version:  2009.04.25
% Create:   2009.04.23
% Coder:    Xin Liang

error(nargchk(4, 5, nargin, 'struct'));
error(nargoutchk(0, 1, nargout, 'struct'));
if nargin == 4 || length(flag) < 2
    flag = bodfltchk(flag, [0.95 0.05], @check);
end
Judge = 0;

p1 = Data.F - Data0.F;
p2 = StepSize * dot(Data0.g, Step);
q1 = dot(Data.g, Step);
q2 = dot(Data0.g, Step);
if  q1 > flag(1) * q2 && p1 <= flag(2) * p2
    Judge = 1;
end

function chk = check(flag)
chk = ones(size(flag));
lf = length(flag);
if lf < 2 || ~isnumeric(flag(2)) || flag(2) >= .5 || flag(2) <= 0
    chk(2) = 0;
end
if lf < 1 || ~isnumeric(flag(1)) || flag(1) >= 1 || flag(1) <= flag(2)
    chk(1) = 0;
end
chk = chk(1:lf);
