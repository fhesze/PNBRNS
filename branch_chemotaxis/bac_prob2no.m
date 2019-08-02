function [bac_mov_no]  = bac_prob2no( bacno,bac_prob )
% Determine the amount of microbes goes to each pore
%   Detailed explanation goes here


% bacno=1000;
% bac_prob=[1 -1 8 0 6 8 9]

bac_prob(bac_prob <0)=0;
Pmax=sum(bac_prob);

rand_no=Pmax.*rand(bacno,1);

prob_ranges=zeros(size(bac_prob));
bac_mov_no=zeros(size(bac_prob));

prob_ranges(1)=bac_prob(1);

for i=1:length(prob_ranges)-1
    prob_ranges(i+1)=prob_ranges(i)+bac_prob(i+1);
    ind_bac=find(rand_no>prob_ranges(i)& rand_no<prob_ranges(i+1));
    bac_mov_no(i+1)=length(ind_bac);
end

bac_mov_no(1)=bacno-sum(bac_mov_no(2:end));


end

