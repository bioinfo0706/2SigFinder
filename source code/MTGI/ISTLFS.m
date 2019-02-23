function [label_sequence ScoringList_ttest_temp]=ISTLFS(Seq,G_signatures,sig,interation_whole,feature_size,each_window_length)
% To reveal potentially dramatic changes at a smaller scale, 
% we propose an iteration of small-scale t-tests with large-scale feature selection (IST-LFS) 
% to quantify the compositional differences of a region from the host
% user¡¯s setting. 
%   Seq-> Input genome
%   G_signatures->Genomic signatures
%   Sig-> Standard deviation of the mean of the window scores to select windows whose scores are large enough to
%       be considered statistically significant.
%   Interation_whole-> Periods of time that are repeated to select windows whose scores are large enough to 
%       be considered statistically significant.
%   feature_size-> Size of selected features by the proposed kurtosis
%   each_window_length-> Length of the windows
%   MTGIpick enables robust identification of genomic islands from a single genome
%   Qi Dai, 20 Apri 2015 

    % combinate larger windows
    mer=cell(1,1);
    mer{1}=G_signatures;
    sig_len=length(sig);
    label_whole_kernal=zeros(sig_len,interation_whole,size(mer{1},1));
    label_whole_ttest=zeros(sig_len,interation_whole,size(mer{1},1));
    % combinate larger windows
    [row column]=size(mer{1});
    for width_window=4:4
        % compute frequency of larger windows
        temp_mer_width=zeros(row-width_window+1,column);
        if width_window==1
            temp_mer_width=mer{1};
        else
           for order_feature=1:size(mer{1},2)
                slide_value=zeros(width_window,width_window+size(mer{1},1)-1);
                for zz=1:width_window
                    slide_value(zz,:)=[zeros(1,width_window-zz) mer{1}(:,order_feature)' zeros(1,zz-1)];
                end
                slide_value(:,1:(width_window-1))=[];
                slide_value(:,(size(mer{1},1)-(width_window-1)+1):size(mer{1},1))=[];
                temp_mer_width(:,order_feature)=mean(slide_value)';
           end
        end
    end
    
 for sig_order=1:sig_len
     fprintf(['Compute whole score with Significant: ' num2str(sig(sig_order)) '\n' ]);
    % interation of process  
    tem_mer_interation_ttest=[];
    tem_mer_interation_ttest=temp_mer_width;
    label_window_ttest=1:size(mer{1},1);
    GI_predict_Interation=zeros(1,size(mer{1},1));
    Label_window_whole=zeros(2,size(mer{1},1));
    for interation_time=1:interation_whole
       % compute score using kernal density or ttest
       ScoringList_ttest_temp=[];
       ScoringList_ttest_temp=WholeComputeScore(tem_mer_interation_ttest,'eandm',feature_size);
       
       % find the significant windows
       position_window_whole=cell(1,2);
       [position_window_whole{2} island_position_window_whole_ttest]=GIwindowposition(ScoringList_ttest_temp,sig(sig_order));
       
       % label of predicted window
       for zzz=2:2
           temp_len_p=length(position_window_whole{zzz});
           for zzzz=1:temp_len_p
               Label_window_whole(zzz,label_window_ttest(position_window_whole{zzz}(zzzz):position_window_whole{zzz}(zzzz)+width_window-1))=1;
           end
       end
       label_whole_ttest(sig_order,interation_time,:)=Label_window_whole(2,:);
       
      % update the frequency of windows
       tem_mer_interation_ttest(position_window_whole{2},:)=[];
       label_window_ttest(position_window_whole{2})=[];
       fprintf(['Iteration: ' num2str(interation_time) '\n' ]);
    end % interation
 end % significant

    % filter the begin and end of predicted label and detele some GI with
    % length less than 7
    for sig_order=1:sig_len
        for interation_time=1:interation_whole
            temp_predict_ttest=reshape(label_whole_ttest(sig_order,interation_time,:),1,size(mer{1},1));
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
                    label_whole_ttest(sig_order,interation_time,island_position_window_ttest(zz,1):island_position_window_ttest(zz,2))=0;
                elseif island_position_window_ttest(zz,2)==size(mer{1},1)
                    label_whole_ttest(sig_order,interation_time,island_position_window_ttest(zz,1):island_position_window_ttest(zz,2))=0;
                elseif island_position_window_ttest(zz,2)-island_position_window_ttest(zz,1)<7
                    label_whole_ttest(sig_order,interation_time,island_position_window_ttest(zz,1):island_position_window_ttest(zz,2))=0;
                end
            end
        end % interation_time
    end % sig_order
          
   % set label to each bases of sequence
   % find star and end position of GI region
     label_sequence=zeros(1,length(Seq));
     tem_position=cumsum(each_window_length);
     tem_label=reshape(label_whole_ttest(sig_len,interation_whole,:),1,size(mer{1},1));
     num_window_total=length( tem_position);
     label_sequence(1,1:tem_position(1,1))=tem_label(1);
     for zz=2:num_window_total
         label_sequence(1,tem_position(1,zz-1)+1:tem_position(1,zz))=tem_label(zz);
     end

