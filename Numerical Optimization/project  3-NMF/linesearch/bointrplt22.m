function [NewSize New] = ...
    bointrplt22(ObjFun, Point, Step, StepSize, Data, varargin)
%BOINTRPLT22 Get stepsize by 2-Point Quadric Intropolation 
%
% In inexact line search problem P: argmin ff(a) = f(x+a*d), along with g
% the gradient of f, when just getting StepSize a0, this function can get a
% better one. Actually it is the minimal point of 
%           Q(a) = q0*a^2 + q1*a + q2
% which satisfy terms below:
%   Q(0) = ff(0) = f(x), Q'(0) = ff'(0) = g'*d, Q(a0) = ff(a0) = f(x+a0*d)
%
% Call
% [NewSize New] =  bointrplt22(ObjFun, Point, Step, StepSize, Data)
% [NewSize New] =  bointrplt22(ObjFun, Point, Step, StepSize, Data,
%                                                       varargin)
%
% Input & Output
% just look up the help of M-function bolinesearch 

% Version:  2009.05.10
% Create:   2009.04.24
% Coder:    Xin Liang

error(nargchk(5, inf, nargin, 'struct'));
error(nargoutchk(0, 2, nargout, 'struct'));

[Q g] = feval(ObjFun, Point, varargin{:});
Qp = dot(g, Step);
Qa = Data.F;

NewSize = - .5 * Qp * StepSize^2 / ( Qa - Q - Qp*StepSize );
if NewSize < 0
    NewSize = StepSize;
end
[New.F New.g] = feval(ObjFun, Point + NewSize*Step, varargin{:});  
