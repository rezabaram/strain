clear all

% plots fitness or N distribution for each time step

init_file=30000; %6570;
num_files=44133; %6570+330; 
step=10;

s=1;

for j=init_file:step:num_files % number of lines in the file "fitness"

filename=['line',int2str(j)]

a=load(filename);

%time(j)=a(1);

most_fit(s)=max( a(2:length(a)) );
least_fit(s)=min( a(2:length(a)) );

s=s+1;

end

max_fitness=max(most_fit);
min_fitness=min(least_fit);

bins=min_fitness:(max_fitness-min_fitness)/100:max_fitness;
s=1;

for j=init_file:step:num_files % number of lines in the file "fitness"

filename=['line',int2str(j)]

a=load(filename);

%time(j)=a(1);

h=figure(1)
hist( a(2:length(a)) , bins)
axis ([-10 10 0 40])
axis 'auto x'
set(gca,'nextplot','replacechildren');
title(sprintf('file name %s, fitness distribution for t=%g years',filename,a(1)))
xlabel('fitness values')
ylabel('counts')

%saveas(h,['fitness_distribution_snapshot',int2str(s),'.fig'])
%Mov(s)=getframe(gcf); 
s=s+1;


%title(sprintf('R_0 = %g, 1/nu = %g days, i(0) = %g',R0,1/nu*365,i0))

%pause(0.25)

end

%movie2avi(Mov, 'fitness_distribution.avi')

%saveas(h,['fitness_distribution_snapshot','.eps']);
%saveas(h,['fitness_distribution_snapshot','.fig']);