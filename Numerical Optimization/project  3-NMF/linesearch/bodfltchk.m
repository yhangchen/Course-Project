function optcrct = bodfltchk(opt, default, specialrule)
%BODLTCHK Correct the unfit opt to default
%
% This function can only correct the Inf/NaN automatically.
%
% Call
% optcrct = bodfltchk(opt, default)
% optcrct = bodfltchk(opt, default, specialrule)
% 
% Input
% specialrule:  vector or function_handle
%                   The function should be declared as a M-function
%                       chk = specialrule(opt)
%                   and chk have the same length as opt.
%                   So do the vector.

% Version:  2009.04.24
% Create:   2009.04.24
% Coder:    Xin Liang

error(nargchk(2, 3, nargin, 'struct'));
error(nargoutchk(0, 1, nargout, 'struct'));
optcrct = default;

chk = isfinite(opt);
if nargin == 3
    if isa(specialrule,'function_handle')
        chk = chk & specialrule(opt);
    else
        chk = chk & specialrule;
    end
end

for i = 1:min(length(opt),length(default))
    if chk(i)
        optcrct(i) = opt(i);
    end
end
