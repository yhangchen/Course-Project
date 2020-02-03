function [U, V, P] = DGS_2(U0, V0, P0, F, G, A, B1, B2, D, N, v1)
% DGS in the slides "V_cycle"
	U=U0;V=V0;P=P0;h = 1/N;
    G_S = tril(A);
    for iter = 1:v1
        U = U + G_S\(F - A*U - B1*P);
        V = V + G_S\(G - A*V - B2*P);
        u = reshape(U,N-1,N);
        v = reshape(V,N-1,N)';
        p = reshape(P,N,N);
        d = reshape(D,N,N);
        for i=2 : N-1
            for j=2 : N-1
                r = - d(i,j) -(u(i,j)-u(i-1,j)+v(i,j)-v(i,j-1))/h;
                delta = r * h/4;
                u(i-1,j) = u(i-1,j) - delta;
                u(i,j) = u(i,j) + delta;
                v(i,j-1) = v(i,j-1) - delta;
                v(i,j) = v(i,j) + delta;
                p(i,j) = p(i,j) + r;
                p(i+1,j) = p(i+1,j) - r/4;
                p(i-1,j) = p(i-1,j) - r/4;
                p(i,j+1) = p(i,j+1) - r/4;
                p(i,j-1) = p(i,j-1) - r/4;
            end
        end

        i = 1;
        for j=2:N-1
            r = - d(i,j) - (u(i,j)+v(i,j)-v(i,j-1))/h;
            delta = r * h / 3;
            v(i,j-1) = v(i,j-1) - delta;
            u(i,j) = u(i,j) + delta;
            v(i,j) = v(i,j) + delta;
            p(i,j) = p(i,j) + r;
            p(i+1,j) = p(i+1,j) - r/3;
            p(i,j+1) = p(i,j+1) - r/3;
            p(i,j-1) = p(i,j-1) - r/3;
        end

        j = 1;
        for i=2:N-1
            r = -d(i,j)-(u(i,j)-u(i-1,j)+v(i,j))/h;
            delta = r*h/3;
            u(i-1,j) = u(i-1,j) - delta;
            u(i,j) = u(i,j) + delta;
            v(i,j) = v(i,j) + delta;
            p(i,j) = p(i,j) + r;
            p(i+1,j) = p(i+1,j) - r/3;
            p(i-1,j) = p(i-1,j) - r/3;
            p(i,j+1) = p(i,j+1) - r/3;
        end
        
        i = N;
        for j=2:N-1
            r = -d(i,j) -(-u(i-1,j)+v(i,j)-v(i,j-1))/h;
            delta = r*h/3;
            u(i-1,j) = u(i-1,j) - delta;
            v(i,j-1) = v(i,j-1) - delta;
            v(i,j) = v(i,j) + delta;
            p(i,j) = p(i,j) + r;
            p(i-1,j) = p(i-1,j) - r/3;
            p(i,j-1) = p(i,j-1) - r/3;
            p(i,j+1) = p(i,j+1) - r/3;
        end

        j = N;
        for i=2:N-1
            r = -d(i,j)-(u(i,j)-u(i-1,j)-v(i,j-1))/h;
            delta = r*h/3;
            u(i-1,j) = u(i-1,j) - delta;
            u(i,j) = u(i,j) + delta;
            v(i,j-1) = v(i,j-1) - delta;
            p(i,j) = p(i,j) + r;
            p(i+1,j) = p(i+1,j) - r/3;
            p(i-1,j) = p(i-1,j) - r/3;
            p(i,j-1) = p(i,j-1) - r/3;
        end
        
        %(1,1)
        r = -d(1,1)-(u(1,1)+v(1,1))/h;
        delta = r*h/2;
        u(1,1) = u(1,1)+delta;
        v(1,1) = v(1,1)+delta;
        p(1,1) = p(1,1)+r;
        p(2,1) = p(2,1)-r/2;
        p(1,2) = p(1,2)-r/2;
        
        %(1,N)
        r = -d(1,N)-(u(1,N)-v(1,N-1))/h;
        delta = r*h/2;
        u(1,N) = u(1,N)+delta;
        v(1,N-1) = v(1,N-1)-delta;
        p(1,N) = p(1,N)+r;
        p(2,N) = p(2,N)-r/2;
        p(1,N-1) = p(1,N-1)-r/2;
        
        %(N,1)
        r = -d(N,1)-(-u(N-1,1)+v(N,1))/h;
        delta = r*h/2;
        u(N-1,1) = u(N-1,1)-delta;
        v(N,1) = v(N,1)+delta;
        p(N,1) = p(N,1)+r;
        p(N-1,1) = p(N-1,1)-r/2;
        p(N,2) = p(N,2)-r/2;
        
        %(N,N)
        r = -d(N,N)-(-u(N-1,N)-v(N,N-1))/h;
        delta = r*h/2;
        u(N-1,N) = u(N-1,N)+delta;
        v(N,N-1) = v(N,N-1)+delta;
        p(N,N) = p(N,N)+r;
        p(N-1,N) = p(N-1,N)-r/2;
        p(N,N-1) = p(N,N-1)-r/2;


        U = u(:);
        v_T = v';
        V = v_T(:);
        P = p(:);
    end
end