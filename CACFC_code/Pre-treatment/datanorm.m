function data=datanorm(data,method)
if strcmp(method,'range')%min-max
    data_min = min(data,[],1);
    data_max = max(data,[],1);
    data = (data-data_min)./(data_max-data_min);
elseif strcmp(method,'var')%z-score(x-mean)/std
    data = (data-mean(data))./sqrt(std(data).^2+1e-6);
    %data = (data-mean(data))./std(data);
elseif strcmp(method,'kar')%
    data = data./sum(data,1);
elseif strcmp(method,'mean')%
    data = (data-mean(data));
 elseif strcmp(method,'norm')%
     [nSmp,~] = size(data);
     for i = 1:nSmp
         data(i,:) = data(i,:) ./ max(1e-12,norm(data(i,:)));
     end
end
data(isnan(data))=0;
end