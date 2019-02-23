function [GI_sequence_left_position]=Semtgi(Seq,G_signature,mer_moment,label_window,feature_size,sig,Iteration,d_feature,sliding_step,each_window_length)
% Delete the predicted genomic islands and update all of windows of the
% genome; then, predict genomic island again using the MTGIpick algorithm
% user's setting. 
%   label_window-> Results from the ILST-DSFS
%   Seq-> Input genome
%   G_signature-> Genome sigatures
%   mer_moment-> Higher moments of the genomic signatures
%   feature_size-> The size of selected genomic signatures
%   Sig-> Standard deviation of the mean of the window scores to select windows whose scores are large enough to
%       be considered statistically significant.
%   Iteration-> Periods of time that are repeated to select windows whose scores are large enough to 
%       be considered statistically significant.
%   d_feature-> The size of selected dynamic genomic signatures
%   sliding_step-> The following continued windows to create the ith large sliding window

%   MTGIpick enables robust identification of genomic islands from a single genome
%   Qi Dai, 20 Apri 2015 

%% detect GI from the rest region
 mer{1}=G_signature;
 left_window_mer=mer_moment{2}(label_window,:);
% compute score using kernal density or ttest
 ScoringList_ttest_temp=[];
 ScoringList_ttest_temp=WholeComputeScore(left_window_mer,'eandm',feature_size);
% find the significant windows
 position_window_left=cell(1,2);
 [position_window_left{2} island_position_window_left_ttest]=GIwindowposition(ScoringList_ttest_temp,sig);
 Label_window_left=zeros(1,size(mer{1},1));
 Label_window_left(1,label_window(position_window_left{2}))=1; 
% filter the begin and end of predicted label and detele some GI with
  temp_predict_ttest= Label_window_left;
 
% find the star and end position window with value 1
  %ttest
  tem_window_isindex_com=find(temp_predict_ttest==1);
  island_position_window_ttest=zeros(1,2);
  tt=1;
  segment_label=1;
  while segment_label<=length(tem_window_isindex_com) 
        extend_window=0;
        while segment_label+extend_window<length(tem_window_isindex_com)&tem_window_isindex_com(segment_label+1+extend_window)-tem_window_isindex_com(segment_label+extend_window)<2
              extend_window=extend_window+1;
        end
        island_position_window_ttest(tt,1)=tem_window_isindex_com(segment_label);
        island_position_window_ttest(tt,2)=tem_window_isindex_com(segment_label+extend_window);
        segment_label=segment_label+extend_window+1;
        tt=tt+1;
   end
   for zz=1:tt-1
      if island_position_window_ttest(zz,1)==1
         Label_window_left(1,island_position_window_ttest(zz,1):island_position_window_ttest(zz,2))=0;
      elseif island_position_window_ttest(zz,2)==size(mer{1},1)
         Label_window_left(1,island_position_window_ttest(zz,1):island_position_window_ttest(zz,2))=0;
      end
   end
          
   % set label to each bases of sequence
   % find star and end position of GI region
     label_sequence_left=zeros(1,length(Seq));
     tem_position=cumsum(each_window_length);
     tem_label=Label_window_left;
     num_window_total=length( tem_position);
     label_sequence_left(1,1:tem_position(1,1))=tem_label(1);
     for zz=2:num_window_total
         label_sequence_left(1,tem_position(1,zz-1)+1:tem_position(1,zz))=tem_label(zz);
     end
     
    % interation of process based on the first dynamic feature 
    tem_mer_interation_left=cell(1,5);
    for zz=1:5
        tem_mer_interation_left{zz}=mer_moment{zz}(label_window,:);
    end
    label_window_left=1:size(mer{1},1);
    GI_sequence_left_position=zeros(1,2);
    GI_index_left=1;
for interation_time=1:Iteration
    % select the main window
    % compute the score using the first 20 dynamic features
    fprintf(['Left Iteration: ' num2str(interation_time) '\n']);
    ScoringList_kernal_dynamic=zeros(5,size(tem_mer_interation_left{1},1));
    for order_moment=1:5    
        ScoringList_kernal_dynamic(order_moment,:)=DynamicComputeScore(tem_mer_interation_left{order_moment},'kernal',sliding_step,20);
    end
    
    % find GI windows     
    temp_window_isindex_kernal=cell(1,6);
    temp_window_isindex_ttest=cell(1,6);
    for order_moment=1:5 
       [temp_window_isindex_kernal{order_moment} temp_a]=GIwindowposition(ScoringList_kernal_dynamic(order_moment,:),0.05);
   end
    temp_window_isindex_kernal{6}=unique([temp_window_isindex_kernal{1} temp_window_isindex_kernal{2} temp_window_isindex_kernal{3} temp_window_isindex_kernal{4} temp_window_isindex_kernal{5}]);
   
     % find the star and end position window with value 1
     tem_window_isindex_com=label_window(label_window_left(temp_window_isindex_kernal{6}));
     if length(tem_window_isindex_com)>0
     island_position_window=zeros(1,2);
     tt=1;
     segment_label=1;
     while segment_label<=length(tem_window_isindex_com) 
          extend_window=0;
          while segment_label+extend_window<length(tem_window_isindex_com)&tem_window_isindex_com(segment_label+1+extend_window)-tem_window_isindex_com(segment_label+extend_window)<2
               extend_window=extend_window+1;
          end
          island_position_window(tt,1)=tem_window_isindex_com(segment_label);
          island_position_window(tt,2)=tem_window_isindex_com(segment_label+extend_window);
          segment_label=segment_label+extend_window+1;
          tt=tt+1;
     end

     % find star and end position of GI region
     tem_position=cumsum(each_window_length);
     window_star=[];window_end=[];
     window_star=[1 tem_position(1:length(each_window_length)-1)+1];
     window_end=[tem_position(1:length(each_window_length))];
     num_GI_predicted=size(island_position_window,1);
     island_position=zeros(num_GI_predicted,2);
     for order_GI_predicted=1:num_GI_predicted
        island_position(order_GI_predicted,1)=window_star(island_position_window(order_GI_predicted,1));
        island_position(order_GI_predicted,2)=window_end(island_position_window(order_GI_predicted,2));
     end

     %  segmentation of sequence using CG content
     num_region_GI=size(island_position_window,1);
     for order_region_GI=1:num_region_GI
         temp_sequence_GI=Seq(island_position(order_region_GI,1):island_position(order_region_GI,2));
          % segment sequence using MJSD
          test_seq1=temp_sequence_GI;
          score_entropy=[];
          %score_entropy=ntentropy(test_seq1,'window',1000);
          score_entropy=ntentropy(test_seq1);
          negative=find(score_entropy<0);
          positive=find(score_entropy>0);
          index_position=zeros(1,length(test_seq1));
          index_position(negative)=-1;
          index_position(positive)=1;
          % use value to index the change position
          index_change=index_position(1:length(test_seq1)-1).*index_position(2:length(test_seq1));
          index_position=find(index_change==-1);
          % find the start and end positions of segmentation
          tem_change=[1 index_position length(temp_sequence_GI)];
          total_number_change=length(tem_change);
          segment_position=zeros(1,2);
          tt=1;
          segment_label=1;
          while segment_label<length(tem_change) 
                extend_window=1;
                while segment_label+extend_window<length(tem_change)&tem_change(segment_label+extend_window)-tem_change(segment_label+extend_window-1)<1000
                      extend_window=extend_window+1;
                end
          segment_position(tt,1)=tem_change(segment_label);
          segment_position(tt,2)=tem_change(segment_label+extend_window);
          segment_label=segment_label+extend_window;
          tt=tt+1;
          end
          % compute the avearge score for each segmentation
          segment_sequence=segment_position+island_position(order_region_GI,1)-1;
          tem_score=zeros(1,size(segment_sequence,1));
          for zz=1:size(segment_sequence,1)
              tem_score(1,zz)=sum(label_sequence_left(1,segment_sequence(zz,1):segment_sequence(zz,2)-1))/length(label_sequence_left(1,segment_sequence(zz,1):segment_sequence(zz,2)-1));
          end
          value_index=find(tem_score>0.15);
          if length(value_index)==1
             island_sequence_position=zeros(1,2);
             tt=1;
             island_sequence_position(tt,1)=value_index;
             island_sequence_position(tt,2)=value_index;
          else
          island_sequence_position=zeros(1,2);
          tt=1;
          segment_label=1;
             while segment_label<length(value_index) 
                extend_window=0;
                while segment_label+extend_window<length(value_index)&value_index(segment_label+extend_window+1)-value_index(segment_label+extend_window)<2
                      extend_window=extend_window+1;
                end
                island_sequence_position(tt,1)=value_index(segment_label);
                island_sequence_position(tt,2)=value_index(segment_label+extend_window);
                segment_label=segment_label+extend_window+1;
                tt=tt+1;
             end
          end
          if length(value_index)>0
              for zz=1:size(island_sequence_position,1)
                  GI_sequence_left_position(GI_index_left,1)=segment_sequence(island_sequence_position(zz,1),1);
                  GI_sequence_left_position(GI_index_left,2)=segment_sequence(island_sequence_position(zz,2),2);
                  GI_index_left=GI_index_left+1;
              end
          end  
        end
     else
         break;
     end
 
     % update the main window and feature selection
     for order_moment=1:5
         tem_mer_interation_left{order_moment}(temp_window_isindex_kernal{6},:)=[];
     end
     label_window_left(temp_window_isindex_kernal{6})=[];
end