clear all

% plots fitness or N distribution for each time step

init_file=1;
num_files=16647;
nbins=50;

most_fit=zeros(1,num_files);
least_fit=zeros(1,num_files);
n=zeros(1,nbins+1);

for j=init_file:1:num_files % number of lines in the file "fitness"

filename=['linef',int2str(j)]

a=load(filename);

fitness_differences=a(3:end)-a(2);

most_fit(j)=max( fitness_differences );
least_fit(j)=min( fitness_differences );

end

max_fitness=max(most_fit);
min_fitness=min(least_fit);

bins=min_fitness:(max_fitness-min_fitness)/nbins:max_fitness;

for j=init_file:1:num_files % number of lines in the file "fitness"

filename=['linef',int2str(j)]

a=load(filename);

fitness_differences=a(3:end)-a(2);

%h=figure(1);
%axes('FontSize',16);
n=n+hist(fitness_differences,bins);
%axis([min_fitness max_fitness 0 125]);
%title('distribution of fitness differences','fontsize',16);
%set(gca,'nextplot','replacechildren');
%xlabel('values of fitness differences','fontsize',16);
%ylabel('counts','fontsize',16);

%saveas(h,['fitness_distribution_snapshot',int2str(j-init_file+1),'.fig'])

%Mov(j-init_file+1)=getframe(); 

%title(sprintf('R_0 = %g, 1/nu = %g days, i(0) = %g',R0,1/nu*365,i0))

%pause(0.25)

end

n=n/num_files;
h=figure(1)
axes('FontSize',16);
bar(bins,n)
title('time averaged distribution of f_i-<f>','fontsize',16);
xlabel('f_i-<f>','fontsize',16);
ylabel('counts','fontsize',16);

saveas(h,['time_averaged_distribution_of_fitness_differences.eps']);
saveas(h,['time_averaged_distribution_of_fitness_differences.fig']);