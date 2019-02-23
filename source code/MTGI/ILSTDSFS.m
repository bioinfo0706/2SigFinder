function [predict_label GI_sequence_whole_position label_window]=ILSTDSFS(Seq,mer_moment,feature_size,sliding_step,interation_time,label_sequence,each_window_length)
% At a large scale, we investigate the variability of higher moments of each tetranucleotide and 
% designe an iteration of large-scale statistical testing using dynamic signals from small-scale feature 
% selection (ILST-DSFS), to identify large, multi-window segments.
% user¡¯s setting. 
%   Seq-> Input genome
%   mer_moment-> Higher moments of the genomic signatures
%   feature_size-> The size of selected dynamic genomic signatures
%   Sig-> Standard deviation of the mean of the window scores to select windows whose scores are large enough to
%       be considered statistically significant.
%   Iteration_time-> Periods of time that are repeated to select windows whose scores are large enough to 
%       be considered statistically significant.
%   sliding_step-> The following continued windows to create the ith large sliding window
%   label_sequence->  The results from the IST-LFS method
%   each_window_length-> Length of the windows

%   MTGIpick enables robust identification of genomic islands from a single genome
%   Qi Dai, 20 Apri 2015 


% interation of process based on the first dynamic feature 
    tem_mer_interation=cell(1,5);
    tem_mer_interation=mer_moment;
    predict_label=zeros(size(mer_moment{1},1),1);
    label_window=1:size(mer_moment{1},1);
    GI_predict_Interation_kernal=zeros(1,size(mer_moment{1},1));
    GI_sequence_whole_position=zeros(1,2);
    GI_index_whole=1;
for interation_time=1:interation_time
    % select the main window
    % compute the score using the first 20 dynamic features
    fprintf(['Iteration: ' num2str(interation_time) '\n']);
    ScoringList_kernal_dynamic=zeros(5,size(tem_mer_interation{1},1));
    for order_moment=1:5    
        ScoringList_kernal_dynamic(order_moment,:)=DynamicComputeScore(tem_mer_interation{order_moment},'kernal',sliding_step,feature_size);
    end
    
    % find GI windows     
    temp_window_isindex_kernal=cell(1,6);
    for order_moment=1:5 
       [temp_window_isindex_kernal{order_moment} temp_a]=GIwindowposition(ScoringList_kernal_dynamic(order_moment,:),0.05);
   end
    temp_window_isindex_kernal{6}=unique([temp_window_isindex_kernal{1} temp_window_isindex_kernal{2} temp_window_isindex_kernal{3} temp_window_isindex_kernal{4} temp_window_isindex_kernal{5}]);
   
     % find the star and end position window with value 1
     tem_window_isindex_com=label_window(temp_window_isindex_kernal{6});
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
         temp_window_mer=mer_moment{1}(island_position_window(order_region_GI,1):island_position_window(order_region_GI,2),:);
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
              tem_score(1,zz)=sum(label_sequence(1,segment_sequence(zz,1):segment_sequence(zz,2)-1))/length(label_sequence(1,segment_sequence(zz,1):segment_sequence(zz,2)-1));
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
                  GI_sequence_whole_position(GI_index_whole,1)=segment_sequence(island_sequence_position(zz,1),1);
                  GI_sequence_whole_position(GI_index_whole,2)=segment_sequence(island_sequence_position(zz,2),2);
                  GI_index_whole=GI_index_whole+1;
              end
          end  
        end
     else
         break;
     end

      % region associated with GI
       predict_label(label_window(temp_window_isindex_kernal{6}))=[1]; 
       
     % update the main window and feature selection
     for order_moment=1:5
         tem_mer_interation{order_moment}(temp_window_isindex_kernal{6},:)=[];
     end
     label_window(temp_window_isindex_kernal{6})=[];
end