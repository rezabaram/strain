clear all


% plots fitness for each strain

num_files=419;

for j=1:1:num_files %419 % number of lines in the file "fitness"

filename=['line',int2str(j)];

a=load(filename);

%{
figure(1)
plot( a(2:end) , '--rs','LineWidth',2, 'MarkerEdgeColor','k', 'MarkerFaceColor','g', 'MarkerSize',10)
%}

time(j)=a(1);

s=1;

clear nonzero_fitness

for i=2:1:length(a)
    if(a(i)~=0)    
        nonzero_fitness(s)=a(i);
        s=s+1;
    end
end

fittest(j)=max( nonzero_fitness );

end

figure
plot(time,fittest, '--rs','LineWidth',1, 'MarkerEdgeColor','k', 'MarkerFaceColor','g', 'MarkerSize',2 )

xlabel('t (years)')
ylabel('maximum fitness')

title('\chi(d)=1 \forall d')