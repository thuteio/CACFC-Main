function [idx] = CACFC(X,c,nc,lambda,alpha)
N=size(X{1},1);
V=length(X);
inc=randperm(N,nc);

Z=zeros([nc N V]);
J=zeros([nc N V]);
sigma=zeros([nc N V]);

W=zeros([nc N V]);

A=cell([1 V]);
E=cell([1 V]);
Y=cell([1 V]);
D=zeros([1 V]);
gamma=zeros([1 V]);

gold=(sqrt(5)-1)/2;

for i=1:V
    D(i)=size(X{i},2);
    X{i}=X{i}';
    A{i}=X{i}(:,inc);
    A{i}=normalize(A{i},'norm');
    A{i}(isnan(A{i}))=0;
    Y{i}=zeros(size(X{i}));
    E{i}=zeros(D(i),N);
    disM=A{i}'*X{i};
    gamma(i)=mean(pdist2(X{i}',A{i}','squaredeuclidean'),"all");
    W(:,:,i)=softmax(disM/(gamma(i)));
end

sigmai=1e-5;sigmaimax=1e10;
y=1e-5;ymax=1e10;
rou=2;
th=1e-6;

stepi=ones([1 V]);
options = optimoptions('quadprog', 'Algorithm', 'interior-point-convex', 'Display', 'off','StepTolerance',1e-3,'MaxIterations',100);
Aeq=ones([1 nc]);

LA1=@(Xi,Ai,Zi,Ei,Yi,Wi,gammai)(y/2)*norm((Xi-Ai*Zi-Ei+(Yi/y)),'fro')^2-lambda*trace(Zi'*log(Wi));
Lsum=@(Xi,Ai,Zi,Ei,Wi)norm((Xi-Ai*Zi-Ei),'fro')^2;

for iter=1:1e3
    %updata Z
    for i=1:V
        Xb=X{i}-E{i}+Y{i}/y;
        Jb=J(:,:,i)-sigma(:,:,i)/sigmai;
        H=(y/2)*A{i}'*A{i}+sigmai/2*eye(nc);
        F=-(y*Xb'*A{i}+lambda*log(W(:,:,i)')+sigmai*Jb');
        if N<500
            for j=1:N
                if nc<50
                    Z(:,j,i)=QPAS(H, F(j,:)')';
                else
                    Z(:,j,i)=quadprog((H+H'), F(j,:), [], [], Aeq, 1, zeros([nc 1]), [],[],options);
                end
            end
        else
            parfor j=1:N
                if nc<50
                    Z(:,j,i)=QPAS(H, F(j,:)')';
                else
                    Z(:,j,i)=quadprog((H+H'), F(j,:), [], [], Aeq, 1, zeros([nc 1]), [],[],options);
                end
            end
        end
    end

    %updata A
    for i=1:V

        L1=(-y)*(X{i}-A{i}*Z(:,:,i)-E{i}+Y{i}/y)*Z(:,:,i)';
        L2=(lambda/gamma(i))*X{i}*(Z(:,:,i)-W(:,:,i))';
        L=(L1-L2);
        JAold=LA1(X{i},A{i},Z(:,:,i),E{i},Y{i},W(:,:,i));

        for sxr=1:30
            Ai=A{i}-L*stepi(i);
            Ai=normalize(Ai,'norm');
            Ai(isnan(Ai))=0;
            Wi=softmax((Ai'*X{i})/gamma(i));
            JAnew=LA1(X{i},Ai,Z(:,:,i),E{i},Y{i},Wi);
            if JAold>JAnew
                break;
            else
                stepi(i)=stepi(i)*gold;
            end
        end
        if sxr<30
            A{i}=Ai;
            W(:,:,i)=Wi;
        end
    end

    %updata E
    tempE=[];
    for i=1:V
        tempE=[tempE;X{i}-A{i}*Z(:,:,i)+Y{i}/y];
    end
    Econcat=prox21(tempE,2/y);
    ro_b=0;
    E{1}=Econcat(1:size(X{1},1),:);
    ro_end=size(X{1},1);
    for i=2:V
        ro_b=ro_b + size(X{i-1},1);
        ro_end=ro_end + size(X{i},1);
        E{i}=Econcat(ro_b+1:ro_end,:);
    end

    %updata J
    OM=Z+sigma/sigmai;
    J=prox_FDI(OM,alpha/sigmai,c);

    Lnew=0;
    for i=1:V
        Lnew=Lnew+Lsum(X{i},A{i},Z(:,:,i),E{i},W(:,:,i));
    end
    fprintf("iter=%d\tcost=%.6f\n",iter,Lnew);
    if iter>30
        if abs(Lnew-Lold)<th
            break;
        end
    end
    Lold=Lnew;

    %updata Lagrange
    sigma=sigma+sigmai*(Z-J);
    for i=1:V
        Y{i}=Y{i}+y*(X{i}-A{i}*Z(:,:,i)-E{i});
    end
    sigmai=min(rou*sigmai,sigmaimax);
    y=min(rou*y,ymax);
end

Zsum=[];Wsum=[];
for i=1:V
    Zsum=[Zsum Z(:,:,i)'];
    Wsum=[Wsum W(:,:,i)'];
end
[Ui,~,~]=svd((Zsum+Wsum)/2,'econ');
[U1,~,~]=svd((Zsum),'econ');
[U2,~,~]=svd((Wsum),'econ');
[U3,~,~]=svd([Zsum,Wsum],'econ');

idx=[kmeans(Ui(:,1:c),c,'Replicates',25) kmeans(U1(:,1:c),c,'Replicates',25) kmeans(U2(:,1:c),c,'Replicates',25) kmeans(U3(:,1:c),c,'Replicates',25)];
end


