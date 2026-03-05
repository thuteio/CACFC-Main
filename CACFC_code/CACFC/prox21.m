function E = prox21(D,alpha)
dil=size(D,2);
E=zeros(size(D));
for i=1:dil
    dl=D(:,i);
    normdl=norm(dl,2);
    if normdl>alpha
        E(:,i)=((normdl-alpha)/normdl)*dl;
    end
end
end