clc;
clear;
warning off;
addpath('CVI',"Pre-treatment","CACFC")


%filename:webkb; wiki; 3sources; youtube; prokaryotic; reuters3;
filename = "webkb";
X=importdata("data/"+filename+".mat");
data=X.data;
label=X.label;
c=length(unique(label));
N=size(data{1},1);
data = Normdata(data);

switch filename
    case "webkb"
        nc=1*c;	lambda=1e-4;	alpha=1e-1;
    case "wiki"
        nc=1*c;	lambda=1e-4;	alpha=1e-5;
    case "3sources"
        nc=3*c;	lambda=1e-1;	alpha=1e-3;
    case "youtube"
        nc=6*c;	lambda=1e-6;	alpha=1;
    case "prokaryotic"
        nc=2*c;	lambda=1e-1;	alpha=1e-5;
    case "reuters3"
        nc=1*c;	lambda=1e-5;	alpha=1e-3;
end


idx = CACFC(data,c,nc,lambda,alpha);

[ACC,NMI,ARI,purity]=CVI(idx,label);

fprintf("=====================\nACC = %.4f\nNMI = %.4f\nARI = %.4f\npurity = %.4f\n=====================",ACC,NMI,ARI,purity);
