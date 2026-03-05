function data = Normdata(data)
V=length(data);
N=size(data{1},1);
minD=1e10;
for v=1:V
    minD=min(minD,size(data{v},2));
end
if V>=6
    str="var";
elseif N<=200
    str="norm";
elseif minD<=5
    str="range";
else
    str="no";
end
for v=1:V
    data{v}=datanorm(data{v},str);
end
end