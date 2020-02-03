function V_cycle_Uzawa_test
% N = [6];
% V = [2,4,8];
% L = [1,2,4];
% Alpha=[0.5,1,1.5,2];

% N = [7,8,9];
% V = [2,3,4];
% L = [1,2,4];
% Alpha = [0.5,1,1.5]; 
% 
% N = [10,11];
% V = [2,3];
% L = [1,2];
% Alpha = [1]; 


for i = 1:length(N)
for k1 = 1:length(Alpha)
for k = 1:length(L)
for j = 1:length(V)
if L(k)<N(i)
choice=[];
choice.alpha = Alpha(k1);
choice.smoother = 'exact';
choice.max_iter = 1000;
choice.L = L(k);
choice.v1 = V(j);
choice.v2 = V(j);
choice.eps = 1e-8;
choice.n = N(i);
out = V_cycle(choice);
fprintf('%d,%.1f, %d, %d, %d, %.2f, %.2e, %.2e, %.2e\n',2^N(i),Alpha(k1), L(k), V(j), out.iter, out.time, out.res, out.u_err, out.p_err)
end
end
end
end
end
end


