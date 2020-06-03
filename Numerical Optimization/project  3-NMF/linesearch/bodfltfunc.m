function opt = bodfltfunc(opt, default)
%BODLTCHK Fill default function_handle to opt when blank
%
% Call
% optcrct = bodfltfunc(opt, default)

% Version:  2009.04.25
% Create:   2009.04.24
% Coder:    Xin Liang


error(nargchk(2, 2, nargin, 'struct'));
error(nargoutchk(0, 1, nargout, 'struct'));

if ~isa(opt, 'function_handle')
    opt = default;
end
