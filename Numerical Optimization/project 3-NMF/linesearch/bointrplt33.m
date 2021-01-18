function [NewSize New] = ...
    bointrplt33(ObjFun, Point, Step, StepSize, Data, varargin)
%BOINTRPLT33 Get stepsize by 3-Point Cubic Intropolation 
%
% In inexact line search problem P: argmin ff(a) = f(x+a*d), along with g
% the gradient of f, when just getting StepSize [a1 a0], this function can
% get a better one. Actually it is the minimal point of 
%           Q(a) = q0*a^3 + q1*a^2 + q2*a + q3
% which satisfy terms below:
%       Q(0) = ff(0) = f(x), Q'(0) = ff'(0) = g'*d
%       Q(a0) = ff(a0) = f(x+a0*d), Q(a1) = ff(a1) = f(x+a1*d)
% Notice that if StepSize is a scalar or there is only a0 without a1, this
% function will call bointrplt22 automatically.
%
% Call
% [NewSize New] =  bointrplt33(ObjFun, Point, Step, StepSize, Data)
% [NewSize New] =  bointrplt33(ObjFun, Point, Step, StepSize, Data,
%                                                       varargin)
%
% Input & Output
% just look up the help of M-function bolinesearch 

% Version:  2009.05.17
% Create:   2009.04.24
% Coder:    Xin Liang

error(nargchk(5, inf, nargin, 'struct'));
error(nargoutchk(0, 2, nargout, 'struct'));

if length(StepSize) == 1
    [NewSize New] = ...
        bointrplt22(ObjFun, Point, Step, StepSize, Data, varargin{:});
    NewSize = [NewSize StepSize];
    New = [New Data];
    return;
end

[Q g] = feval(ObjFun, Point, varargin{:});
Qp = dot(g, Step);
Qa0 = Data(2).F;
Qa1 = Data(1).F;
a0 = StepSize(2);
a1 = StepSize(1);

ff = 0;
if a0 ~= 0 && a1 ~= 0 && a1 ~= a0
    p = [a0^2 -a1^2; -a0^3 a1^3] * [Qa1-Q-Qp*a1; Qa0-Q-Qp*a0] ...
        / ( a0^2 * a1^2 *(a1-a0) );
    p = roots([p.*[3; 2]; Qp]);
    if isreal(p)
        p = max(p);
        NewSize = [p StepSize(1)];
        ff = 1;
    end
end
if ff ~= 1
    [NewSize New] = ...
        bointrplt22(ObjFun, Point, Step, StepSize(1), Data, varargin{:});
    NewSize = [NewSize StepSize(1)];
end
[New.F New.g] = feval(ObjFun, Point + NewSize(1)*Step, varargin{:});
New = [New Data(1)];

    
