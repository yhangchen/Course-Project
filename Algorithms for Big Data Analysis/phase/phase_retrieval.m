function test_err_1 = phase_retrieval(n, choice, L, flag)

% modified from https://viterbi-web.usc.edu/~soltanol/WFcode.html, combine
% two parts of codes together, and we modified the operator for CDP.

if nargin < 4
    flag = 1;
end
if nargin < 3
    L = 0;
end
if strcmp(choice,'gaussian')
    %% Make 1D Gaussian signal and data 
    x = randn(n,1) + 1i*randn(n,1);
    if L == 0
        L = 4.5;
    end
    m = round(L*n);                     
    A0 = 1/sqrt(2)*randn(m,n) + 1i/sqrt(2)*randn(m,n);
    y = abs(A0*x).^2;
    A = @(I) A0*I;  % Input is n x 1 signal, output is m x 1 array
    At = @(Y) A0'*Y/m;            % Input is m x 1 array, output is n x 1 signal          

elseif strcmp(choice, 'cdp')
    %% Make signal
    x = randn(n,1) + 1i*randn(n,1);
    %% Make masks and linear sampling operators
    if L == 0
        L = 6;% Number of masks
    else
        if fix(L) ~= L
            error('L should be integer');
        end
    end
    % Sample phases: each symbol in alphabet {1, -1, i , -i} has equal prob. 
    Masks = randsample([1i -1i 1 -1],n*L,true);
    Masks = reshape(Masks,n,L);
    % Sample magnitudes and make masks 
    temp = rand(size(Masks));
    Masks = Masks .* ( (temp <= 0.2)*sqrt(3) + (temp > 0.2)/sqrt(2) );
    % Make linear operators; A is forward map and At its scaled adjoint (At(Y)*numel(Y) is the adjoint) 
    A = @(I)  vec(fft(conj(Masks) .* repmat(I,[1 L])));  % Input is n x 1 signal, output is nL x 1 array
    At = @(Y) mean(Masks .* ifft(reshape(Y,size(Masks))), 2);            % Input is nL x 1 array, output is n x 1 signal          
    % Data 
    y = abs(A(x)).^2;  
else
    error('Wrong type')
end
%% test
max_iter = 2500;
tau0 = 330;
tic
[z, Err] = Wirtinger(A, At, y, x, max_iter, tau0);
test_err_1 = Err(end)
t1 = toc

% when n is large, comment the following part.
% tic
% [z2, err] = PhaseLift_CVX(A, At, y, x);
% test_err_2 = err
% t2 = toc
if flag == 1
    % plot Wirtinger flow.
    figure
    semilogy(Err,'linewidth',1)
    title(string(choice)+'  L = '+string(L))
    set(gca,'FontSize',12);
    xlabel('iteration')
    ylabel('relative error')
    grid on
end

end