function [StepSize info perf] = ...
    bolinesearch(ObjFun, Point, Step, Rule, varargin)
%BOLINESEARCH Find the answer to problem P: a = argmin f(x+a*d)
%
% This function will use exact or inexact line search to solve the problem.
% When doing inexact line search, you can choose Armijo-Goldstein, Wolfe,
% strong Wolfe or you own criterion.
% Record g as the gradient of f. 
%
% Call
% [StepSize err] = bolinesearch(ObjFun, Point, Step)
% [StepSize err] = bolinesearch(ObjFun, Point, Step, Rule)
% [StepSize err] = bolinesearch(ObjFun, Point, Step, Rule, p1,p2,...)
% [StepSize err perf] = bolinesearch(......)
%
% Input 
% ObjFun:   f & g in P, function_handle
%               The function should be declared as a M-function
%                   [F g] = f(x, p1,p2,...)
%               F is a scalar, along with g an n-vector.
% Point:    x in P, n-vector
% Step:     d in P, n-vector
% Rule:     option & method & criterion to solve P, struct
%       Rule.crtr:      criterion, function_handle
%               The function should be declared as a M-function
%                   Judge = criterion(Step, StepSize, Data0, Data, flag)
%               Judge is a logical number with 1 perfect, and Step is just
%               as metioned above StepSize below, and Data0 & Data are
%               structs including possible fields F g, flag is a array
%               including several parameters in the criterion.
%           choice: boarmgld, bowlf, bostwlf
%       Rule.mthd:      method to get new point, function_handle
%               The function shoule be declared as a M-function
%                   [NewSize New] = method(ObjFun, Point, Step, StepSize,  
%                                               Data, p1,p2,...)
%               notice here StepSize & NewSize can be a scalar array(to
%               make coding easy), but only the 1st element is just needed
%               actually. And in that situation the Data & New will be
%               struct arrays.
%           choice: bointrplt22, bointrplt33
%       Rule.opt:       options of iteration, scalar array
%           opt(1):     0 - exact line search, use .618 method
%                       else - inexact, need crtr & mthd
%              (2): upper bound of a
%              (3): maximum of iterations
%            (4:5): criterion flag, but also can be .618 method flag which
%                   have one element e with default 1e-3
%       default:        bostwlf - bointrplt33 - [1 10 10 0.95 0.05]
%                                       or [0 10 25 1e-4] (if opt(1) = 0)
%
% Output
% StepSize: a in P just the answer, scalar
% info:     status to execute this function, integral array
%       info(1):    exit code
%               0 - Successful call
%               1 - Reach the maximum of iterations
%              -1 - a is not real
%           (2):    number of iterations
%           (3):    number of evaluating ObjFun
% perf:     other useful data, struct
%       perf.x:     new point after iteration
%       perf.F:     function value at new point
%       perf.g:     gradient at new point

% Version:  2009.05.19
% Create:   2009.04.22
% Coder:    Xin Liang

error(nargchk(3, inf, nargin, 'struct'));
error(nargoutchk(0, 3, nargout, 'struct'));
if nargin == 3 || isempty(Rule)
    Rule.opt = 1;
    Rule.crtr = -1;
    Rule.mthd = -1;
end
if ~Rule.opt(1)
    Rule.opt = bodfltchk(Rule.opt, [0 10 25 1e-4]);
else
    Rule.opt = bodfltchk(Rule.opt, [1 10 10 0.95 0.05]);
    Rule.crtr = bodfltfunc(Rule.crtr, @bostwlf);
    Rule.mthd = bodfltfunc(Rule.mthd, @bointrplt33);
end
if isempty(Step)
    Step = zeros(size(Point));
end
StepSize = Rule.opt(2);
info = [0 0 0];

if ~Rule.opt(1)
    k = ( sqrt(5) - 1 ) / 2;
    lb = 0;
    ub = Rule.opt(2);
    fl = 0;
    fr = 0;
    while ub - lb >= Rule.opt(4)*Rule.opt(2)
        if ~fl
            al = lb + (1-k) * (ub - lb);
            Fl = feval(ObjFun, Point + al * Step, varargin{:});
            info(3) = info(3) + 1;
        end
        if ~fr
            ar = lb + k * (ub - lb);
            Fr = feval(ObjFun, Point + ar * Step, varargin{:});
            info(3) = info(3) + 1;
        end
        if Fl < Fr
            ub = ar;
            ar = al;
            Fr = Fl;
            fr = 1;
            fl = 0;
        else
            lb = al;
            al = ar;
            Fl = Fr;
            fr = 0;
            fl = 1;
        end
        info(2) = info(2) + 1;
        if info(2) >= Rule.opt(3)
            info(1) = 1;
            break;
        end
    end
    StepSize = (lb + ub) / 2;
    perf.x = Point + StepSize*Step;
    [perf.F perf.g] = feval(ObjFun, perf.x, varargin{:});
    [F1 g1] = feval(ObjFun, Point, varargin{:});
    if F1 < perf.F
        perf.x = Point;
        perf.F = F1;
        perf.g = g1;
        StepSize = 0;
    end
    [F1 g1] = feval(ObjFun, Point + Rule.opt(2)*Step, varargin{:});
    if F1 < perf.F
        perf.x = Point;
        perf.F = F1;
        perf.g = g1;
        StepSize = Rule.opt(2);
    end
    info(3) = info(3) + 3;
else
   [Data0.F Data0.g] = feval(ObjFun, Point, varargin{:});
   [Data.F Data.g] = feval(ObjFun, Point + Rule.opt(2)*Step, varargin{:});
   info(3) = info(3) + 2;
   Judge = feval(Rule.crtr, Step, StepSize, Data0, Data, Rule.opt(4:5));
   while ~Judge
       [StepSize Data] = feval(Rule.mthd, ...
           ObjFun, Point, Step, StepSize, Data, varargin{:});
       info(2) = info(2) + 1;
       info(3) = info(3) + 1;
       if ~isreal(StepSize(1)) 
           info(1) = -1;
           break;
       end
       if StepSize(1) < 0
           StepSize(1) = 0;
       end
       Judge = feval(Rule.crtr, ...
           Step, StepSize(1), Data0, Data(1), Rule.opt(4:5));
       if info(2) >= Rule.opt(3)
           info(1) = 1;
           break;
       end
   end
   StepSize = StepSize(1);
   perf.x = Point + StepSize*Step;
   perf.F = Data(1).F;
   perf.g = Data(1).g;
end
