function Plot_DGS
N = [6,7,8,9,10,11];
V = [4];
NmL = [1];


for i = 1:length(N)
for k = 1:length(NmL)
for j = 1:length(V)
if NmL(k)<N(i)
choice=[];
choice.smoother = '2';
choice.max_iter = 1000;
choice.L = N(i)-NmL(k);
choice.v1 = V(j);
choice.v2 = V(j);
choice.eps = 1e-8;
choice.n = N(i);
out = V_cycle(choice);
u = out.U;
v = out.V;
p = out.P;
n = 2^N(i);
diff = [u;v]-mat_U(n);
diffp = p - mat_P(n);
diff1 = reshape(diff(1:n*(n-1)),n-1,n);
diff2 = reshape(diff(n*(n-1)+1:2*n*(n-1)),n-1,n);
diff2 = diff2';
diff3 = reshape(diffp,n,n);
h_fig1 = figure('Visible', 'on');
s=pcolor(diff1);
s.FaceColor='interp';
set(s,'LineStyle','none');
title(['N=',num2str(n),' error of U']);
zdir=[0,0,1];
xlabel('1');
ylabel('2');
colorbar;
saveas(h_fig1, [num2str(n),'U','.png']);
h_fig1 = figure('Visible', 'off');

h_fig1 = figure('Visible', 'on');
s=pcolor(diff2);
s.FaceColor='interp';
set(s,'LineStyle','none');
title(['N=',num2str(n),' error of V']);
zdir=[0,0,1];
xlabel('1');
ylabel('2');
colorbar;
saveas(h_fig1, [num2str(n),'V','.png']);
h_fig1 = figure('Visible', 'off');          

h_fig1 = figure('Visible', 'on');
s=pcolor(diff3);
s.FaceColor='interp';
set(s,'LineStyle','none');
title(['N=',num2str(n),' error of P']);
zdir=[0,0,1];
xlabel('1');
ylabel('2');
colorbar;
saveas(h_fig1, [num2str(n),'P','.png']);
h_fig1 = figure('Visible', 'off');               

end
end
end
end
end
