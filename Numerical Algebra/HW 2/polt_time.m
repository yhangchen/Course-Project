x = zeros(10,1);
t = zeros(10,1);
err = zeros(10,1);
for N = 1:10
    x(N)=10*N-1;
    [u1,t(N)]=poisson_choice(@f,10*N-1,'gauss');
    err(N) = norm(u1-real_sol(10*N-1),inf);
end
hold off
plot(x,t)%这里我们仅仅给出绘制时间代码，绘制误差代码可改为
%plot(x,err)
hold on
scatter(x,t)
xlabel('每一方向网格数目')
ylabel('时间/秒')
title('Gauss消去')%随需要而改变。
hold off
