function J = prox_FDI(OM, lambda,c)
[n1,n2,n3]=size(OM);
n12=min(n1,n2);
Yf = fft(OM, [], 3);

Yf(isnan(Yf)) = 0;
Yf(isinf(Yf)) = 0;

Uf = zeros(n1, n12, n3);
Vf = zeros(n2, n12, n3);
Sf = zeros(n12,n12, n3);

h=c*log(c);

for i=1:n3
    [Uf(:,:,i), Sf(:,:,i), Vf(:,:,i)] = svd(Yf(:,:,i), 'econ');
    s = diag(Sf(:, :, i));
    for o=1:length(s)
        sold=s(o);
        a=(o-c)/c;
        for k=1:1e3
            snew=max(0,(s(o)-lambda*((a*h^(a^2))/(a+sold).^2) ));
            if(snew-sold)<1e-3
                break;
            end
            sold=snew;
        end
        s(o)=snew;
        if snew==0
            s(o:end)=0;
            break;
        end
    end
    Sf(:, :, i) = diag(s);
    Yf(:,:,i)=Uf(:,:,i)*Sf(:, :, i)*Vf(:, :, i)';
end

J=ifft(Yf,[],3);

end