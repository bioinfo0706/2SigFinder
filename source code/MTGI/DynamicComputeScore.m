
function  ScoringList=DynamicComputeScore(mer_whole,method,sliding_step,feature_size)
   
     ScoringList=[];
    %% select the main window
    core_window_whole=[];
    variance_window_whole=[];
    variance_window_whole=var(mer_whole')';
    densitiesSig1_var_whole=[];
    densitiesSig2_var_whole=[];
    [densitiesSig1_var_whole,densitiesSig2_var_whole,xmaxlist,ymaxlist]=CalculateDensities(variance_window_whole, 1,0);% '5sigma' or '3sigma'    
    ScoringList_var_whole=[];
    ScoringList_var_whole=ComputeScore(densitiesSig1_var_whole,densitiesSig2_var_whole,variance_window_whole);
    core_window_whole=mer_whole;
    core_window_whole(find(ScoringList_var_whole==1),:)=[];
    
    %% Estimate kernal density function 
     if strmatch(method,'kernal')
       %% compute the score using kernal density function
        densitiesSig1_whole=[];
        densitiesSig2_whole=[];
        [densitiesSig1_whole,densitiesSig2_whole,xmaxlist,ymaxlist]=CalculateDensities(core_window_whole, 3,4);% '5sigma' or '3sigma'\
    end
    
    %% Select features along the sliding windows using kertosis and skewness
    length_section=size(mer_whole,1);
    ScoringList=zeros(1,length_section);
    sliding_window_position=1:sliding_step:length_section;
    dynamic_feature=zeros(length(sliding_window_position),2,feature_size);
    for num_region_select=1:length(sliding_window_position)
            if num_region_select>length(sliding_window_position)-1
               tem_mer_region_select=[];
               sta_tem_mer=sliding_window_position(num_region_select);
               end_tem_mer=length_section;
            else
               tem_mer_region_select=[];
               sta_tem_mer=sliding_window_position(num_region_select);
               end_tem_mer=sliding_window_position(num_region_select+1)-1;
            end

            %selct the feature
            tem_mer_region_select=mer_whole(sta_tem_mer:end_tem_mer,:); %dynamic region
            num_feature=size(mer_whole,2);
            ks_value=zeros(2,num_feature);
            for i1=1:num_feature
                ks_value(1,i1)= kurtosis(tem_mer_region_select(:,i1)); % 
                ks_value(2,i1) = skewness(tem_mer_region_select(:,i1));  % 
            end
            ks_ascending=zeros(2,num_feature);
            ks_position=zeros(2,num_feature);
            for i1=1:2
                [ks_ascending(i1,:) ks_position(i1,:)]=sort(ks_value(i1,:),'descend');
            end
            dynamic_feature(num_region_select,1,:)=ks_position(1,1:feature_size);
            dynamic_feature(num_region_select,2,:)=ks_position(2,1:feature_size);
            
            % compute score
            if strmatch(method,'kernal')
            % compute the score using kernal density function
            ScoringList_kernal_local=[];
            ScoringList_kernal_local=ComputeScore1(densitiesSig1_whole,densitiesSig2_whole,tem_mer_region_select,ks_position(1,1:feature_size));%largest
            ScoringList(1,sta_tem_mer:end_tem_mer)=ScoringList_kernal_local;
            else
            % compute the score using eye and main t-test method
            ScoringList_ttest_local=[];
            ScoringList_ttest_local=eandmvalue1(tem_mer_region_select,core_window_whole,5,ks_position(1,1:feature_size));
            ScoringList(1,sta_tem_mer:end_tem_mer)=ScoringList_ttest_local;
          end
    end
     
  