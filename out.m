%n=2; % prints N
n=3; % prints fitness

Ncutoff=40;

for j=1:1:9

filename=['single0000',int2str(j)];

a=load(filename);

if (length(a)~=0)
    if (max(a(:,2))>Ncutoff)
        figure(n)
        j
        plot(a(:,1),a(:,5)-a(:,6))
        %plot(a(:,1),a(:,n)-a(:,4))
        %plot(a(:,1),a(:,n))
        %plot(a(:,1),a(:,4))
        hold on
    end
end
end


for j=10:1:99

filename=['single000',int2str(j)];

a=load(filename);

if (length(a)~=0) 
    if (max(a(:,2))>Ncutoff)
    figure(n)
    j
    plot(a(:,1),a(:,5)-a(:,6))
    %plot(a(:,1),a(:,n)-a(:,4))
    %plot(a(:,1),a(:,n))
    %plot(a(:,1),a(:,4))
    hold on
    end
end

end

for j=100:1:999

filename=['single00',int2str(j)];

a=load(filename);

if (length(a)~=0) 
    if (max(a(:,2))>Ncutoff)
    figure(n)
    j
    plot(a(:,1),a(:,5)-a(:,6))
    %plot(a(:,1),a(:,n)-a(:,4))
    %plot(a(:,1),a(:,n))
    %plot(a(:,1),a(:,4))
    hold on
    end
end

end
