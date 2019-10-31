x = zeros(91,1);
t = zeros(91,1);
err = zeros(91,1);
for N = 1:91
    x(N)=N+8;
    [u1,t(N)]=poisson_choice(@f,N+8,'gauss_band');
    err(N) = norm(u1-real_sol(N+8),inf);
end
hold off
plot(x,t)%这里我们仅仅给出绘制时间代码，绘制误差代码可改为
%plot(x,err)
hold on
scatter(x,t)
xlabel('每一方向网格数目')
ylabel('时间/秒')
title('带状Gauss消去')%随需要而改变。
hold off