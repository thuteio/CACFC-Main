function ARI_value = ARI(true_labels, cluster_labels)


n = length(true_labels);

a = 0;
b = 0;
c = 0;
d = 0;

for i = 1 : n
    for j = i+1 : n
        if (( true_labels(i,1) ==  true_labels(j,1) ) && ( cluster_labels(i,1) ==  cluster_labels(j,1) ))
            a = a+1;
        end
        
        if (( true_labels(i,1) ==  true_labels(j,1) ) && ( cluster_labels(i,1) ~=  cluster_labels(j,1) ))
            b = b+1;
        end
        
        if (( true_labels(i,1) ~=  true_labels(j,1) ) && ( cluster_labels(i,1) ==  cluster_labels(j,1) ))
            c = c+1;
        end
        
        if (( true_labels(i,1) ~=  true_labels(j,1) ) && ( cluster_labels(i,1) ~=  cluster_labels(j,1) ))
            d = d+1;
        end
        
    end %for j = i+1 : n
end %for i = 1 : n

Up = a - (a+b)*(a+c)/(a+b+c+d);

Down = ((a+b)+(a+c))/2 - (a+b)*(a+c)/(a+b+c+d);

ARI_value = Up / Down;