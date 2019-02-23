function   ScoringList=ComputeScore1(densityWMM,densitiesSig2,data,ChosenTetra) %densityWMM=densitiesSig1
tem_score=zeros(size(data,1),1);
if nargin<4
    ChosenTetra=1:size(data,2);
end
Num_feature=length(ChosenTetra);
for l=1:Num_feature 
    % compute score using kernal density function
     CI=FindCI(densityWMM{ChosenTetra(l)},densitiesSig2{ChosenTetra(l)},0.05);%[0.05, 0.025, 0.01]
     for AlphaNumber=1:length(CI)
         sz1=find(data(:,ChosenTetra(l))<CI{AlphaNumber}(1));
         sz2=find(data(:,ChosenTetra(l))>CI{AlphaNumber}(2));
         ssz=zeros(size(data,1),1);
         ssz(sz1)=[1];
         ssz(sz2)=[1];
         tem_score=tem_score+ssz;
     end
end
ScoringList=sum(sum(tem_score));


