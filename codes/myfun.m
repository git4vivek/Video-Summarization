function [keyframe] = myfun(c1)
%PCA
x = c1;
OutputSize = 6;
[Rowsx, Columnsx] = size(x);  % find size of input matrix
m=mean(x);                  % find mean of input matrix
y=x-ones(size(x,1),1)*m;    % normalise by subtracting mean
c=cov(y);                   % find covariance matrix
[V,D]=eig(c);               % find eigenvectors (V) and eigenvalues (D) of covariance matrix
[D,idx] = sort(diag(D));    % sort eigenvalues in descending order by first diagonalising eigenvalue matrix, idx stores order to use when ordering eigenvectors
D = transpose(D(end:-1:1));
V = V(:,idx(end:-1:1));     % put eigenvectors in order to correspond with eigenvalues
V2d=V(:,1:OutputSize);        % (significant Principal Components we use, OutputSize is input variable)
prefinal=transpose(V2d)*transpose(y);
final=transpose(prefinal);            % final is normalised data projected onto eigenspace

%delaunay
deltri = delaunayn(final);   
[Rows_deltri, Columns_deltri] = size(deltri);

%matrix initialised to 0
for i = 1:Rowsx
    for j = 1:Rowsx
        eucl_dist(i,j) = 0;
    end
end

%finding euclidean distances between images obtained from delaunay
for i = 1:Rowsx
    for j = 1:Rows_deltri
        for k = 1:Columns_deltri
            if deltri(j,k) == i
                for l = 1:Columns_deltri
                    if deltri(j,l) ~= i
                        if eucl_dist(i,deltri(j,l)) == 0
                            eucl_dist(i,deltri(j,l)) = (pdist2(final(i),final(deltri(j,l)))*10);
                            eucl_dist(deltri(j,l),i) = eucl_dist(i,deltri(j,l));
                        end
                    end
                end
            end
        end
    end
end


countn=0;
sumn = 0;
sumd = 0;
localdev = 0;

for i = 1:Rowsx
    for j = 1:Rowsx
        if eucl_dist(i,j) ~= 0
            countn= countn+1;
            sumn = sumn + eucl_dist(i,j);
        end
    end
    local_mlength(i) = sumn / countn;
    for k = 1:Rowsx
        if eucl_dist(i,k) ~= 0
            sumd = sumd + (local_mlength(i) - eucl_dist(i,k))*(local_mlength(i) - eucl_dist(i,k));
        end
    end
    localdev(i) = sqrt (sumd/countn);
    sumn =0;
    sumd =0;
    countn = 0;
end

sumn = 0;
for i= 1: Rowsx
    sumn = sumn + localdev(i);
end

global_dev = sumn/Rowsx;

for i = 1:Rowsx
    for j = 1:Rowsx
        if eucl_dist(i,j) ~= 0
            if ((eucl_dist(i,j)) > ( local_mlength(i) - global_dev)) 
                %sprintf('eculidean dist:\n%s', eucl_dist(i,j))
%                 sprintf('threshhold :\n%s', local_mlength(i) - global_dev)
                eucl_dist(i,j) = 0;
                eucl_dist(j,i) = 0; 
            end
        end
    end
end

for i = 1:Rowsx
    for j = 1:Rowsx
        visited(i,j)=0;
    end
end

t=1;

k=0;
cluster = 0;
for i = 1:Rowsx
    flag = 0;
    k = 0;
    for m = 1:Rowsx
        if visited(i,m) == 1
            flag = 1;
        end
    end
    stack = 0;
    if flag == 0   
        cluster(t) = i;
        t = t+1;
        for j = 1:Rowsx
            if (eucl_dist(i,j) ~= 0) &&  (visited(i,j) == 0)
                    k = k+1;
                    stack(k)= j;
                    visited(i,j) = 1;
                    visited(j,i) = 1;
            end
        end
        while k ~= 0
             x = stack(k);
             k=k-1;
             for l = 1:Rowsx
                 if (eucl_dist(x,l) ~= 0) && (visited(x,l)==0)
                         k = k+1;
                         stack(k)= l;
                         visited (x,l)=1;
                         visited (l,x)=1;
                         visited (i,l)=1;
                         visited (l,i)=1;
                  end
             end
        end
    end
end    

for i = 1:Rowsx
    for j = 1:Rowsx
        if eucl_dist(i,j) == 0
            eucl_dist(i,j) = (pdist2(final(i),final(j))*10);
            eucl_dist(j,i) = eucl_dist(i,j);
        end
    end
end
    


k = 1;
index = 0;
[r c] = size(cluster);
while k <= c
        q = cluster(k);
        t = 1;
        clus(t) = q;
        for i = 1:Rowsx
                if visited(q,i) == 1
                    t = t+1;
                    clus(t)=i;
                end     
        end
        [dr cr] = size(clus);
       
        if cr > 1
            for j = 1:cr
                min = 1000000;
                index = 0;
                sum = 0;
                for l= 1:cr
                    sum = sum + eucl_dist(clus(j),clus(l));
                end
                if (sum < min) 
                    min = sum;
                    index = j;
                end
            end
            keyframe(k) = index;
            else
                keyframe(k) = q;
        end
        k = k+1;
        clus=[];
end
clearvars -except keyframe
keyframe = sort(keyframe);
