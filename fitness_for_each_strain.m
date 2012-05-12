clear all

% plots fitness for each strain

num_files=59;

for j=1:1:num_files % number of lines in the file "fitness"

filename=['line',int2str(j)]

a=load(filename);

figure(1)
plot( a(3:end) , '--rs','LineWidth',2, 'MarkerEdgeColor','k', 'MarkerFaceColor','g', 'MarkerSize',10)

end