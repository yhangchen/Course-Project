function b = gauss_band(A,b,m)
%´ø¿í2m+1
n = length(A);
for k = 1:n-1
    if A(k,k) ~=0
        s = min(k+m,n);
        A(k+1:s,k)=A(k+1:s,k)/A(k,k);
        A(k+1:s,k+1:s)=A(k+1:s,k+1:s)-A(k+1:s,k)*A(k,k+1:s);
    end
    if A(k,k)==0
        ME = MException('Singular');
        throw(ME);
    end
end
for j = 1:n-1
    s = min(j+m,n);
    b(j+1:s)=b(j+1:s)-b(j)*A(j+1:s,j);
end
for j = n:-1:2
    s = max(j-m,1);
    b(j)=b(j)/A(j,j);
    b(s:j-1)=b(s:j-1)-b(j)*A(s:j-1,j);
end
b(1)=b(1)/A(1,1);
end