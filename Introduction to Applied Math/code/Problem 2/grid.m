function u = grid(m,n,v,r)
%Solve problem: Au_{n+1} = u_n, where r is the depth and v(1),v(2) specify 
%iteration steps of Gauss-Seidel method.
%Return results of the last layer.
    function U = restriction(le)
        %Restriction matrix from (le-1)x(le-1) to (le/2-1)x(le/2-1)
        U = sparse(zeros((le/2-1)^2,(le-1)^2));
        for i = 1:le/2-1
            for j = 1:le/2-1
                U(((i-1)*(le/2-1)+j),((le-1)*(2*i-2)+(2*j-1)))=1/16;
                U(((i-1)*(le/2-1)+j),((le-1)*(2*i-2)+(2*j+1))) =1/16;
                U(((i-1)*(le/2-1)+j),(le-1)*(2*i)+(2*j-1))=1/16;
                U(((i-1)*(le/2-1)+j),(le-1)*(2*i)+(2*j+1))=1/16;
                U(((i-1)*(le/2-1)+j),(le-1)*(2*i-1)+(2*j-1))=1/8;
                U(((i-1)*(le/2-1)+j),(le-1)*(2*i-1)+(2*j+1))=1/8;
                U(((i-1)*(le/2-1)+j),(le-1)*(2*i-2)+2*j)=1/8;
                U(((i-1)*(le/2-1)+j),(le-1)*(2*i)+2*j)=1/8;
                U(((i-1)*(le/2-1)+j),(le-1)*(2*i-1)+(2*j))= 1/4;
            end
        end
    end

%     function U = improvement(le)
%         U = sparse(zeros((le-1)^2,(le/2-1)^2));
%         %i=j=0
%         U(1,1)=1/4;
%         %i=j=le/2-1
%         U((2*(le/2-1)-1)*(le-1)+2*(le/2-1), ((le/2-1)-1)*(le/2-1)+(le/2-1))=1;
%         U((2*(le/2-1))*(le-1)+2*(le/2-1), ((le/2-1)-1)*(le/2-1)+(le/2-1))=1/2;
%         U((2*(le/2-1)-1)*(le-1)+2*(le/2-1)+1, ((le/2-1)-1)*(le/2-1)+(le/2-1))=1/2;
%         U((2*(le/2-1))*(le-1)+2*(le/2-1)+1, ((le/2-1)-1)*(le/2-1)+(le/2-1))=1/4;
%         %i=0,j=le/2-1
%         U(2*(le/2-1), (le/2-1))=1/2;
%         U(2*(le/2-1)+1, (le/2-1))=1/4;
%         %j=0,i=le/2-1
%         U((2*(le/2-1)-1)*(le-1)+1, ((le/2-1)-1)*(le/2-1)+1)=1/2;
%         U((2*(le/2-1))*(le-1)+1, ((le/2-1)-1)*(le/2-1)+1)=1/4;
%         for j = 1:le/2-2%i=0
%             U(2*j,j) = 1/2;
%             U(2*j+1,j) = 1/4;
%             U(2*j+1,j+1) = 1/4;
%             U((2*(le/2-1)-1)*(le-1)+2*j, ((le/2-1)-1)*(le/2-1)+j)=1;
%             U((2*(le/2-1))*(le-1)+2*j, ((le/2-1)-1)*(le/2-1)+j)=1/2;
%             U((2*(le/2-1)-1)*(le-1)+2*j+1, ((le/2-1)-1)*(le/2-1)+j)=1/2;
%             U((2*(le/2-1)-1)*(le-1)+2*j+1, ((le/2-1)-1)*(le/2-1)+j+1)=1/2;
%             U((2*(le/2-1))*(le-1)+2*j+1, ((le/2-1)-1)*(le/2-1)+j)=1/4;
%             U((2*(le/2-1))*(le-1)+2*j+1, ((le/2-1)-1)*(le/2-1)+j+1)=1/4;
%         end
%         for i = 1:le/2-2%j=0
%             U((2*i-1)*(le-1)+1,(i-1)*(le/2-1)+1)=1/2;
%             U((2*i)*(le-1)+1,(i-1)*(le/2-1)+1)=1/4;
%             U((2*i)*(le-1)+1,i*(le/2-1)+1)=1/4;
%             U((2*i-1)*(le-1)+2*(le/2-1), (i-1)*(le/2-1)+(le/2-1))=1;
%             U((2*i)*(le-1)+2*(le/2-1), (i-1)*(le/2-1)+(le/2-1))=1/2;
%             U((2*i)*(le-1)+2*(le/2-1), (i)*(le/2-1)+(le/2-1))=1/2;
%             U((2*i-1)*(le-1)+2*(le/2-1)+1, (i-1)*(le/2-1)+(le/2-1))=1/2;
%             U((2*i)*(le-1)+2*(le/2-1)+1, (i-1)*(le/2-1)+(le/2-1))=1/4;
%             U((2*i)*(le-1)+2*(le/2-1)+1, (i)*(le/2-1)+(le/2-1))=1/4;
%         end
%         for i = 1:le/2-2
%             for j = 1:le/2-2
%                 U((2*i-1)*(le-1)+2*j, (i-1)*(le/2-1)+j)=1;
%                 U((2*i)*(le-1)+2*j, (i-1)*(le/2-1)+j)=1/2;
%                 U((2*i)*(le-1)+2*j, (i)*(le/2-1)+j)=1/2;
%                 U((2*i-1)*(le-1)+2*j+1, (i-1)*(le/2-1)+j)=1/2;
%                 U((2*i-1)*(le-1)+2*j+1, (i-1)*(le/2-1)+j+1)=1/2;
%                 U((2*i)*(le-1)+2*j+1, (i-1)*(le/2-1)+j)=1/4;
%                 U((2*i)*(le-1)+2*j+1, (i)*(le/2-1)+j)=1/4;
%                 U((2*i)*(le-1)+2*j+1, (i-1)*(le/2-1)+j+1)=1/4;
%                 U((2*i)*(le-1)+2*j+1, (i)*(le/2-1)+j+1)=1/4;
%             end
%         end
% 
%         for j = 1:le/2-2
%             U((le-3)*(le-1)+2*j, (le/2-2)*(le/2-1)+j)=1;
%             U((le-3)*(le-1)+2*j+1, (le/2-2)*(le/2-1)+j)=1;
%         end
%     end

%Gauss-Seidel method with iteration times 'itr'
    function u1 = gs(Ah,uh, Fh,itr)
        Lo = tril(Ah);
        U = -triu(Ah,1);
        u1 = uh;
        for i0 = 1:itr
            u1 = Lo\(U*u1+Fh);
        end
    end

k = 1/m;h = 1/n;
u = zeros((n-1)^2,1);
%Initiate u.
for i = 1:n-1
    for j = 1:n-1
        u((n-1)*(i-1)+j) = sin(pi*i/n)*sin(pi*j/n);
    end
end
k = 1/m;h = 1/n;
u = zeros((n-1)^2,1);
for i = 1:n-1
    for j = 1:n-1
        u((n-1)*(i-1)+j) = sin(pi*i/n)*sin(pi*j/n);
    end
end
%construct A.
for i = 1:(n-1)^2
    dig(i) = 1+4*k/h^2;
end
for i = 1:n^2-2*n
    subdig(i) = -k/h^2;
end
for i = 1:(n-1)*(n-2)
    subsubdig(i) = -k/h^2;
end
x_axis = [1:(n-1)^2,2:(n-1)^2,1:n^2-2*n,1:(n-1)*(n-2),n:(n-1)^2];
y_axis = [1:(n-1)^2,1:n^2-2*n,2:(n-1)^2,n:(n-1)^2,1:(n-1)*(n-2)];
value = [dig,subdig,subdig,subsubdig,subsubdig];
A = sparse(x_axis,y_axis,value);
for i = 1:n-2
    A((n-1)*i+1,(n-1)*i) = 0;
    A((n-1)*i,(n-1)*i+1) = 0;
end

%Store data used repetitively.
%res stores the restriction matrix(sparse)
store = cell(r+1,1);store_A = cell(r+1,1);
store_F = cell(r+1,1);res = cell(r,1);
Ah = A;
store_A{1} = A;
for i1 = 1:r
    res{i1} =  restriction(n/(2^(i1-1)));
end
for i1 = 1:r
    U0 = res{i1};
    Ah = U0*Ah*U0'*4;
    store_A{i1+1} = Ah;
end


tic
count = 0;
for j1 = 1:m
    f = u;
    store_F{1} = f;u0=u;
    init_err = norm(f-A*u0);err = init_err;
    while err > 10^(-6)*init_err
        u0 = gs(A,u0,f,v(1));
        r0 = f-A*u0;
        store{1}=u0;
        for i1 = 1:r
            U0 = res{i1};
            r0 = U0*r0;
            store_F{i1+1} = r0;
            Ah = store_A{i1+1};
            u0 =  gs(Ah,zeros(length(Ah),1),r0,v(1));
            store{i1+1} = u0;
            r0 = -Ah*u0+r0;
        end
        for i1 = r:-1:1
            U0 = res{i1};
            store{i1} = store{i1}+ U0'*u0*4;
            u0 = store{i1};
            u0 =  gs(store_A{i1},store{i1},store_F{i1},v(2));
        end
        err = norm(f-A*u0);
        count = count+1;
    end
    u = u0;
end
toc
X_real = zeros((n-1)^2,1);
for j = 1:(n-1)
    for p = 1:(n-1)
        X_real((n-1)*(j-1)+p) = sin(pi*p/n)*sin(pi*j/n)*exp(-2*pi^2);
    end
end
err = norm(u-X_real)
end
