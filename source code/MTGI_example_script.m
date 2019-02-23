% Qi Dai
% College of Life Sciences, Zhejiang Sci-Tech University, Hangzhou 310018, China
% Department of Biological Sciences, Center for Systems Biology, University of Texas at Dallas,
% Richardson, TX 75080, USA
%
% November 2014
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Warranty Disclaimer and Copyright Notice
% 
% Copyright (C) 2014-2015 Zhejiang Sci-Tech University, Hangzhou 310018, China
% 
% The Zhejiang Sci-Tech University and the authors make no representation about the suitability or accuracy of this software for any purpose, and makes no warranties, either express or implied, including merchantability and fitness for a particular purpose or that the use of this software will not infringe any third party patents, copyrights, trademarks, or other rights. The software is provided "as is". The Institute for Systems Biology and the authors disclaim any liability stemming from the use of this software. This software is provided to enhance knowledge and encourage progress in the scientific community. 
% 
% This is free software; you can redistribute it and/or modify it under the terms of the GNU Lesser General Public License as published by the Free Software Foundation.
% 
% You should have received a copy of the GNU Lesser General Public License
% along with this library; if not, write to the Free Software Foundation,
% College of Life Sciences, Zhejiang Sci-Tech University, Hangzhou 310018, China
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;close all;clc;
addpath(genpath(pwd));
addpath(genpath([pwd '/MTGI/']));

%% Path to fast convolution
   if ispc
       addpath([pwd '\Win_CONVNFFT_Folder']);rmpath([pwd '\Unix_CONVNFFT_Folder']);
   elseif isunix
       addpath([pwd '/Unix_CONVNFFT_Folder']);rmpath([pwd '/Win_CONVNFFT_Folder']);
   end

%% load dataset
    infile=['data/S. Typhi CT18.txt'];    
    fprintf(['Reading genome data' ]);
    fprintf('\n');
    Name=[];Seq=[];
    [Name, Seq] = fastaread(infile);

%% Extract genomic signature
    % Parameters
    k=4;  % word length
    window=1000; % length of the window
    slidelen=1000; % length of the sliding window
    % Calculate the genomic signatures
    mer=[];
    each_window_length=[];
    fprintf(['compute ' num2str(k) '-mer \n']);
    [mer each_window_length]=cmer(Seq,window,slidelen,k);
    
%% Calculate the score of each window using the IST-LFS algorithm
   fprintf(['IST-LFS Algorithm \n' ]);
   % Parameters
   sig=0.05; % Standard deviation of the mean of the window scores to select windows whose scores are large enough to
             % be considered statistically significant
   Iteration=3; % Periods of time that are repeated to select windows whose scores are large enough to 
             % be considered statistically significant.
   feature_size=256; % Size of selected features by the proposed kurtosis
   % Calculate the scores using the IST-LFS algorithm
   label_sequence_s=[];
   Score_ISTLFS=[];
   [label_sequence_s Score_ISTLFS]=ISTLFS(Seq,mer,sig,Iteration,feature_size,each_window_length);
   
%% Identify large segments using ILST-DSFS based on the tetranucleotide frequencies 
%% as well as their normalized first, second, third and fourth standardized moments;
   
   fprintf(['ILST-DSFS Algorithm \n' ]);
   % compute four moments based on a sliding window
    mer_moment=cell(1,5);
    [mer_moment{1} mer_moment{2} mer_moment{3} mer_moment{4} mer_moment{5}]=mermoment(mer,4);
    
   % Parameters
   feature_size=4; % Size of selected dynamic genomic signatures
   sig=0.05; % Standard deviation of the mean of the window scores to select windows whose scores are large enough to
             % be considered statistically significant
   Iteration=10; % Periods of time that are repeated to select windows whose scores are large enough to 
             % be considered statistically significant.
   sliding_step=100; % The following continued windows to create the ith large sliding window
             
  % Identify large, multi-window segments using the ILST-DSFS 
  predict_label=[];label_window=[];
  GI_sequence_whole_position=[];
  [predict_label_l GI_sequence_whole_position label_window]=ILSTDSFS(Seq,mer_moment,feature_size,sliding_step,Iteration,label_sequence_s,each_window_length);
   
%% Predict genomic island using MTGIpick from the rest of the genome
 % Parameters
   feature_size=256; % Size of selected genomic signatures
   sig=0.05; % Standard deviation of the mean of the window scores to select windows whose scores are large enough to
             % be considered statistically significant
   Iteration=1; % Periods of time that are repeated to select windows whose scores are large enough to 
             % be considered statistically significant.
   d_feature=20; % Size of selected dynamic genomic signatures
   sliding_step=50; % The following continued windows to create the ith large sliding window  

   % Predict genomic island from the rest of the genome
   [GI_sequence_left_position]=Semtgi(Seq,mer,mer_moment,label_window,feature_size,sig,Iteration,d_feature,sliding_step,each_window_length);

%% got the final predicted GI region
GI_sequence_predict=[];
if size(GI_sequence_left_position,1)==1
    if GI_sequence_left_position(1,1)~=0 & GI_sequence_left_position(1,2)~=0
       GI_sequence_predict=sort([GI_sequence_whole_position; GI_sequence_left_position]);
    else
       GI_sequence_predict=sort([GI_sequence_whole_position]);
    end
elseif size(GI_sequence_left_position,1)>1
   GI_sequence_predict=sort([GI_sequence_whole_position; GI_sequence_left_position]);
end

%% Refining the boundaries of the predicted genomic islands using the CG content and Markov JS divergence

num_GI_predict=size(GI_sequence_predict,1);
GI_boundary=zeros(num_GI_predict,2);
fprintf(['\n']);
fprintf(['          Predicted Genomic Islands                    ''\n']);
fprintf(['Number          Start            End           Length' '\n']);
for order_GI_predict=1:num_GI_predict
    left_start=max(GI_sequence_predict(order_GI_predict,1)-20000,1);
    left_end=min(GI_sequence_predict(order_GI_predict,2),GI_sequence_predict(order_GI_predict,1)+10000);
    right_start=max(GI_sequence_predict(order_GI_predict,1),GI_sequence_predict(order_GI_predict,2)-10000);
    right_end=min(GI_sequence_predict(order_GI_predict,2)+20000,length(Seq));
    test_seq1=Seq(left_start:left_end);
    test_seq2=Seq(right_start:right_end);
    % the boundary
    GI_boundary(order_GI_predict,1)=CGMJD(test_seq1)+left_start;
    GI_boundary(order_GI_predict,2)=CGMJD(test_seq2)+right_start;
    fprintf([ num2str(order_GI_predict) '              ' num2str(GI_boundary(order_GI_predict,1)) '            ' num2str(GI_boundary(order_GI_predict,2)) '             ' num2str(GI_boundary(order_GI_predict,2)-GI_boundary(order_GI_predict,1)+1) '\n']);
end    
