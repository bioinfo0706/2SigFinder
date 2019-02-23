   function [mer_moment1 mer_moment2 mer_moment3 mer_moment4 mer_moment5]=mermoment(mer,eye)
% Calculate the higher moments of each tetranucleotide 
%   mer-> Input genomic signatures
%   eye-> 2*eye window surrounding the ith window as its neighbourhood 
%   MTGIpick enables robust identification of genomic islands from a single genome
%   Qi Dai, 20 Apri 2015 
    neigh_size=4;
    mer_moment=cell(1,5);
    mer_moment{1}=mer;
    [row column]=size(mer);
    mer_moment{2}=zeros(row,column);
    mer_moment{3}=zeros(row,column);
    mer_moment{4}=zeros(row,column);
    mer_moment{5}=zeros(row,column);
    for order_feature=1:size(mer,2)
         slide_value=zeros(2*eye+1,2*eye+size(mer,1));
         for zz=1:2*eye+1
             slide_value(zz,:)=[zeros(1,2*eye+1-zz) mer(:,order_feature)' zeros(1,zz-1)];
         end
         slide_value(:,1:eye)=[];
         slide_value(:,size(mer,1)+1:eye+size(mer,1))=[];
         mer_moment{2}(:,order_feature)=mean(slide_value)';
         mer_moment{3}(:,order_feature)=var(slide_value)';
         mer_moment{4}(:,order_feature)=skewness(slide_value)';
         mer_moment{5}(:,order_feature)=kurtosis(slide_value)';
    end
    mer_moment1=mer_moment{1};
    mer_moment2=mer_moment{2};
    mer_moment3=mer_moment{3};
    mer_moment4=mer_moment{4};
    mer_moment5=mer_moment{5};
    