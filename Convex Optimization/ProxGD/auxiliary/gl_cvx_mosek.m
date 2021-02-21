function [x, iter, out, penal] = gl_cvx_mosek(x0, A, b, mu, opt)
%GL_CVX_MOSEK �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
% cvx version:
[n, l] = size(x0);
mode = opt{1};
if mode == "TV1D"
    r = 1: n - 1;
    v = ones(1, n - 1);
    D = sparse([r, r], [r, r + 1], [-v, v], n - 1, n);
elseif mode == "TV2D"
    r = 1: l - 1;
    v = ones(1, l - 1);
    D1 = sparse([r, r + 1], [r, r], [-v, v], l, l - 1);
    r = 1: l - 1;
    v = ones(1, l - 1);
    D2 = sparse([r, r], [r, r + 1], [-v, v], n - 1, n);
elseif mode == "ind_nuclear"
    R = opt{2};
end

cvx_solver mosek
cvx_begin
    variable x(n, l);
    if mode == "linf"
        minimize(square_pos(norm(A * x - b, 'fro')) / 2 + ...
            mu * max(max(abs(x))));
    elseif mode == "l12"
        minimize(square_pos(norm(A * x - b, 'fro')) / 2 + ...
            mu * sum(norms(x, 2, 2)));
    elseif mode == "elastic"
        alpha=opt{2};
        minimize(square_pos(norm(A * x - b, 'fro')) / 2 + ...
            mu * (alpha*sum(norms(x, 2, 2))+(1-alpha)/2*square_pos(norm(x,'fro'))));
    elseif mode == "l21"
        minimize(square_pos(norm(A * x - b, 'fro')) / 2 + ...
            mu * sum(norms(x', 2, 2)));
    elseif mode == "TV1D"
        minimize(square_pos(norm(A * x - b, 'fro')) / 2 + ...
            mu * norm(D * x, 1));
    elseif mode == "TV2D"
        minimize(square_pos(norm(A * x - b, 'fro')) / 2 + ...
            mu * sum(norms(D2 * x, 1, 1)) +...
            mu * sum(norms(x * D1, 1, 1)));
    elseif mode == "nuclear"
        minimize(square_pos(norm(A * x - b, 'fro')) / 2 + ...
            mu * norm(x, 2));
    elseif mode == "ind_nuclear"
        minimize(square_pos(norm(A * x - b, 'fro')) / 2);
        norm(A, 2) <= R;
    end
cvx_end
iter = cvx_slvitr;
out = norm(A * x - b, 'fro')^2 / 2;
if mode == "linf"
    penal = max(max(abs(x)));
elseif mode == "l12"
    penal = sum(norms(x, 2, 2));
elseif mode == "elastic"
    penal = alpha*sum(norms(x, 2, 2))+(1-alpha)/2*norm(x,'fro')^2;
elseif mode == "l21"
    penal = sum(norms(x', 2, 2));
elseif mode == "TV1D"
    penal = norm(D * x, 1);
elseif mode == "TV2D"
    penal = sum(norms(D2 * x, 1, 1)) + sum(norms(x * D1, 1, 1));
elseif mode == "nuclear"
    penal = norm(x, 2);
end


end

