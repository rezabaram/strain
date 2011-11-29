clear all

% plots fitness or N distribution for each time step

init_file=4000;
num_files=6644; %79027; % fitness shifted1overd

for j=init_file:1:num_files % number of lines in the file "fitness"

filename=['line',int2str(j)]

a=load(filename);

%time(j)=a(1);

s=1;

clear nonzero_fitness

for i=2:1:length(a)
    if(a(i)~=0)    
        nonzero_fitness(s)=a(i);
        s=s+1;
    end
end

most_fit(j)=max( nonzero_fitness );
least_fit(j)=min( nonzero_fitness );

end

max_fitness=max(most_fit);
min_fitness=min(least_fit);

bins=min_fitness:(max_fitness-min_fitness)/100:max_fitness;

for j=init_file:1:num_files % number of lines in the file "fitness"

filename=['line',int2str(j)]

a=load(filename);

%time(j)=a(1);

s=1;

clear nonzero_fitness

for i=2:1:length(a)
    if(a(i)~=0)    
        nonzero_fitness(s)=a(i);
        s=s+1;
    end
end

figure(1)
hist(nonzero_fitness,bins)
title(sprintf('file name %s, fitness distribution for t=%g',filename,a(1)));

%title(sprintf('R_0 = %g, 1/nu = %g days, i(0) = %g',R0,1/nu*365,i0))

%pause(0.25)

end