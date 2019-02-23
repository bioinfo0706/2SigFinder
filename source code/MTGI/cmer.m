
function [merk_c Nlettercount]=cmer(Seq,window,slidelen,k)
% Calculate k-mer of the input data
% user's setting. 
%   Seq-> Input genome
%   window->Size of sliding window
%   slidelen-> Step of sliding windows.
%   k-> Length of k-mer 


%   MTGIpick enables robust identification of genomic islands from a single genome
%   Qi Dai, 20 Apri 2015 





N=length(Seq);
N_word=4^k;
T_region=1;
N_region= ceil((length(Seq)-window+1)/slidelen);
merk_c=zeros(N_region,N_word);
NU='ATCG';
%the size of kmer

% count the content distribution in the window

% count the content distribution in the window
section=1:slidelen:N-window+1;
Nlettercount=zeros(1,length(section));%each element record the length of each sliding window which from 5000 to 6000
for i=1:length(section)-1
    subsequence=[];subsequence1=[];
    subsequence=Seq(section(i):section(i)+window);
    subsequence1=subsequence;
    subsequence1(subsequence1=='N')='';
    extend_word=1; 
    while length(subsequence1)<window & extend_word<100
          subsequence1=Seq(section(i):section(i)+window+extend_word);
          subsequence1(subsequence1=='N')='';
          extend_word=extend_word+1;
    end
    subsequence=Seq(section(i):section(i)+window+extend_word-2);
    Nlettercount(1,i)=length(subsequence); 
    t=1;
    if k==2
        for k1=1:4
            for k2=1:4
                a=[];b=[];
                a=[NU(k1) NU(k2)];
                b=findstr(subsequence,a);
                merk_c(i,t)=(length(b)+1)/(window+4^k);
                t=t+1;
            end
        end
    elseif k==3
        for k1=1:4
            for k2=1:4
                for k3=1:4
                    a=[];b=[];
                    a=[NU(k1) NU(k2) NU(k3)];
                    b=findstr(subsequence,a);
                    merk_c(i,t)=(length(b)+1)/(window+4^k);
                    t=t+1;
                end
            end
        end
    elseif k==4
        for k1=1:4
            for k2=1:4
                for k3=1:4
                     for k4=1:4
                         a=[];b=[];
                         a=[NU(k1) NU(k2) NU(k3) NU(k4)];
                         b=findstr(subsequence,a);
                         merk_c(i,t)=(length(b)+1)/(window+4^k);
                         t=t+1;
                     end
                end
            end
        end
    elseif k==5
        for k1=1:4
            for k2=1:4
                for k3=1:4
                     for k4=1:4
                         for k5=1:4
                             a=[];b=[];
                             a=[NU(k1) NU(k2) NU(k3) NU(k4) NU(k5)];
                             b=findstr(subsequence,a);
                             merk_c(i,t)=(length(b)+1)/(window+4^k);
                             t=t+1;
                         end
                     end
                end
            end
        end
       elseif k==6
        for k1=1:4
            for k2=1:4
                for k3=1:4
                     for k4=1:4
                         for k5=1:4
                             for k6=1:4
                                 a=[];b=[];
                                 a=[NU(k1) NU(k2) NU(k3) NU(k4) NU(k5) NU(k6)];
                                 b=findstr(subsequence,a);
                                 merk_c(i,t)=(length(b)+1)/(window+4^k);
                                 t=t+1;
                             end
                         end
                     end
                end
            end
        end
       elseif k==7
        for k1=1:4
            for k2=1:4
                for k3=1:4
                     for k4=1:4
                         for k5=1:4
                             for k6=1:4
                                 for k7=1:4
                                    a=[];b=[];
                                    a=[NU(k1) NU(k2) NU(k3) NU(k4) NU(k5) NU(k6) NU(k7)];
                                    b=findstr(subsequence,a);
                                    merk_c(i,t)=(length(b)+1)/(window+4^k);
                                    t=t+1;
                                 end
                             end
                         end
                     end
                end
            end
        end
       elseif k==8
        for k1=1:4
            for k2=1:4
                for k3=1:4
                     for k4=1:4
                         for k5=1:4
                             for k6=1:4
                                 for k7=1:4
                                     for k8=1:4
                                         a=[];b=[];
                                         a=[NU(k1) NU(k2) NU(k3) NU(k4) NU(k5) NU(k6) NU(k7) NU(k8)];
                                         b=findstr(subsequence,a);
                                         merk_c(i,t)=(length(b)+1)/(window+4^k);
                                         t=t+1;
                                     end
                                 end
                             end
                         end
                     end
                end
            end
        end
    end
end


% the last window
    t=1;
    subsequence=[];
    subsequence=Seq(section(length(section)):length(Seq));
    Nlettercount(1,length(section))=length(subsequence);
    if k==2
        for k1=1:4
            for k2=1:4
                a=[];b=[];
                a=[NU(k1) NU(k2)];
                b=findstr(subsequence,a);
                merk_c(length(section),t)=(length(b)+1)/(window+4^k);
                t=t+1;
            end
        end
    elseif k==3
        for k1=1:4
            for k2=1:4
                for k3=1:4
                    a=[];b=[];
                    a=[NU(k1) NU(k2) NU(k3)];
                    b=findstr(subsequence,a);
                    merk_c(length(section),t)=(length(b)+1)/(window+4^k);
                    t=t+1;
                end
            end
        end
    elseif k==4
        for k1=1:4
            for k2=1:4
                for k3=1:4
                     for k4=1:4
                         a=[];b=[];
                         a=[NU(k1) NU(k2) NU(k3) NU(k4)];
                         b=findstr(subsequence,a);
                         merk_c(length(section),t)=(length(b)+1)/(window+4^k);
                         t=t+1;
                     end
                end
            end
        end
    elseif k==5
        for k1=1:4
            for k2=1:4
                for k3=1:4
                     for k4=1:4
                         for k5=1:4
                             a=[];b=[];
                             a=[NU(k1) NU(k2) NU(k3) NU(k4) NU(k5)];
                             b=findstr(subsequence,a);
                             merk_c(length(section),t)=(length(b)+1)/(window+4^k);
                             t=t+1;
                         end
                     end
                end
            end
        end
       elseif k==6
        for k1=1:4
            for k2=1:4
                for k3=1:4
                     for k4=1:4
                         for k5=1:4
                             for k6=1:4
                                 a=[];b=[];
                                 a=[NU(k1) NU(k2) NU(k3) NU(k4) NU(k5) NU(k6)];
                                 b=findstr(subsequence,a);
                                 merk_c(length(section),t)=(length(b)+1)/(window+4^k);
                                 t=t+1;
                             end
                         end
                     end
                end
            end
        end
       elseif k==7
        for k1=1:4
            for k2=1:4
                for k3=1:4
                     for k4=1:4
                         for k5=1:4
                             for k6=1:4
                                 for k7=1:4
                                    a=[];b=[];
                                    a=[NU(k1) NU(k2) NU(k3) NU(k4) NU(k5) NU(k6) NU(k7)];
                                    b=findstr(subsequence,a);
                                    merk_c(length(section),t)=(length(b)+1)/(window+4^k);
                                    t=t+1;
                                 end
                             end
                         end
                     end
                end
            end
        end
       elseif k==8
        for k1=1:4
            for k2=1:4
                for k3=1:4
                     for k4=1:4
                         for k5=1:4
                             for k6=1:4
                                 for k7=1:4
                                     for k8=1:4
                                         a=[];b=[];
                                         a=[NU(k1) NU(k2) NU(k3) NU(k4) NU(k5) NU(k6) NU(k7) NU(k8)];
                                         b=findstr(subsequence,a);
                                         merk_c(length(section),t)=(length(b)+1)/(window+4^k);
                                         t=t+1;
                                     end
                                 end
                             end
                         end
                     end
                end
            end
        end
    end 
  Nlettercount(1,length(section))=length(subsequence);

       