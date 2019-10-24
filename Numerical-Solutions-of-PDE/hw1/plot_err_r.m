err = zeros(21,1);
x = 30:3:90;
for i = 30:3:90
    err(i/3-9)= log(err_estimate_r(@f_r,@g_r,@h_r,i,i));
    %ºŸ…ËM=N°£
    x(i/3-9)= log(i);
end
plot(x,err);
hold on
scatter(x,err);
xlabel('log N')
ylabel('log error')
hold off