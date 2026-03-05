function [ACCi,NMIi,ARIi,purityi]=CVI(idx,label)
len=size(idx,2);
ACCs=zeros([0 len]);
NMIs=zeros([0 len]);
ARIs=zeros([0 len]);
puritys=zeros([0 len]);
for i=1:len
    idxl = bestMap(label,idx(:,i));
    ACCs(i) = ACC(label,idxl);
    NMIs(i) = NMI(label,idx(:,i));
    ARIs(i) = ARI(label,idx(:,i));
    puritys(i) = purity(label,idx(:,i));
end
ACCi=max(ACCs,[],2);
NMIi=max(NMIs,[],2);
ARIi=max(ARIs,[],2);
purityi=max(puritys,[],2);
end