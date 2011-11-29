n=2; % prints N
% n=3; % prints fitness

for j=1:1:9

filename=['single0000',int2str(j)]

a=load(filename);

if (length(a)~=0) 
    figure(n)
    plot(a(:,1),a(:,n))
    hold on
end

end


for j=10:1:99

filename=['single000',int2str(j)]

a=load(filename);

if (length(a)~=0) 
    figure(n)
    plot(a(:,1),a(:,n))
    hold on
end

end

for j=100:1:999

filename=['single00',int2str(j)]

a=load(filename);

if (length(a)~=0) 
    figure(n)
    plot(a(:,1),a(:,n))
    hold on
end

end
