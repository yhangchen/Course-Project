err = zeros(70,1);
x = 30:99;
for i = 30:99
    err(i-29)= log(err_estimate(@f,@g,@h,i));
    x(i-29)= log(i);
end
plot(x,err);
hold on
scatter(x,err);
xlabel('log N')
ylabel('log error')
hold off