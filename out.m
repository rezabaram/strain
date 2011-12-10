%n=2; % prints N
n=3; % prints fitness

Ncutoff=50;

for j=1:1:9

filename=['single0000',int2str(j)]

a=load(filename);

if (length(a)~=0)
    if (max(a(:,2))>Ncutoff)
        figure(n)
        plot(a(:,1),a(:,n))
        hold on
    end
end
end


for j=10:1:99

filename=['single000',int2str(j)]

a=load(filename);

if (length(a)~=0) 
    if (max(a(:,2))>Ncutoff)
    figure(n)
    plot(a(:,1),a(:,n))
    hold on
    end
end

end

for j=100:1:999

filename=['single00',int2str(j)]

a=load(filename);

if (length(a)~=0) 
    if (max(a(:,2))>Ncutoff)
    figure(n)
    plot(a(:,1),a(:,n))
    hold on
    end
end

end
