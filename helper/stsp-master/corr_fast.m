function [d,Ctot] = corr_fast(Tpre,Tpost,Ta,Tb,bin)
if length(Tpre)>1e5 || length(Tpost)>1e5;keyboard;end
% Tlist_pre : presynaptic spiking times
% Tlist_post : postsynaptic spiking times
% Ta : starting point of the histogram before the presynaptic onset
% Tb : ending interval 
% bin : number of bins
% plot : 0/1 plot/don't plot
k=1;
Ctot = [];
Tchunk = 1000;
while (k-1)*Tchunk<max(Tpre)
    Tlist_pre = [];
    Tlist_post = [];
    Tlist_pre = Tpre(Tpre>(k-1)*Tchunk & Tpre<k*Tchunk);
    Tlist_post = Tpost(Tpost>(k-1)*Tchunk & Tpost<k*Tchunk);

    if size(Tlist_post,1)==1
        A = repmat(Tlist_post,length(Tlist_pre),1);
        % A: pre * post, fill with post in row
        
    else
        A = repmat(Tlist_post',length(Tlist_pre),1);
    end

    if size(Tlist_pre,1)==1
        B = repmat(Tlist_pre',1,length(Tlist_post));
        % B: pre * post, fill with pre in col
        
    else
        B = repmat(Tlist_pre,1,length(Tlist_post));
    end

    C = A - B ;
    
%     11 21 31
%     12 22 32
    
%     post1-pre1 post2-pre1 post3-pre1
%     post1-pre2 post2-pre2 post3-pre2

    C(C>Tb | C<Ta | C==0) = NaN;
    C = C(:);
    Ctot = [Ctot C'];
    k = k+1;

end
d = hist(Ctot,bin);