MTGIpick
===
MTGIpick implements multiscale statistical algorithm to predict genomic islands from a single genome. It uses small-scale test with large-scale features to score small region deviating from the host and large-scale statistical test with small-scale features to identify multi-window segments for identification of genomic islands.

There are two main components to create the MTGIpick. The third component is optional.
1. Quantifying the compositional differences of each region from the host using IST-LFS algorithm. 
2. Identifying some multi-window segments using ILST-DSFS algorithm.
3. Refining the boundaries of  the genomic islands using MJSD-CG algorithm.

Below we describe these components using the Matlab MTGI_example_script.m. Note that depending on the size of the data and parameters settings computing the MTGIpick can be demanding in terms of memory and CPU time.


## MTGI_example_script.m

The code consists of four parts:

1. Extract genomic signature
2. Calculate the score of each window using the IST-LFS algorithm
3. Identifying large segments using the ILST-DSFS alogrithm
4. Predicting genomic island using MTGIpick from the rest of the genome
5. Refining the boundaries of the predicted genomic islands using the MJSD-CG algorithm

These four parts are run with an example genomic signal in the main script: msr_example_script.m.
Below we describe these four parts in more detail.


### 1. Extract genomic signature
- create the genomic signatures using a window and a sliding step
 [mer each_window_length]=cmer(Seq,window,slidelen,k);
- mer contains the genomic signatures of each windows, and each_window_length denotes the length of each window.

### 2. Calculate the score of each window using the IST-LFS algorithm
- Quantifying the compositional differences of a region from the host using an iteration of small-scale t-tests with large-scale feature selection (IST-LFS). 
 [label_sequence_s Score_ISTLFS]=ISTLFS(Seq,mer,sig,Iteration,feature_size,each_window_length);
- Predicted values of each bases of the genome using a small-scale t-test, and Score_ISTLFS denotes the score of each window.

### 3. Identifying large segments using the ILST-DSFS alogrithm
- At a large scale, we investigate the variability of higher moments of each tetranucleotide and designe an iteration of large-scale statistical testing using dynamic signals from small-scale feature selection (ILST-DSFS), to identify large, multi-window segments.
 [predict_label_l GI_sequence_whole_position label_window]=ILSTDSFS(Seq,mer_moment,feature_size,sliding_step,Iteration,label_sequence_s,each_window_length);
- Predicted values of each bases of the genome using a large-scale test, GI_sequence_whole_position denotes the start and end positions of the predicted multi-window segments, and label_window are the predicted labels of these windows. 

### 4. Predicting genomic island using MTGIpick from the rest of the genome
- Delete the predicted genomic islands and update all of windows of the genome; then, predict genomic island again using the MTGIpick algorithm
 [GI_sequence_left_position]=Semtgi(Seq,mer,mer_moment,label_window,feature_size,sig,Iteration,d_feature,sliding_step,each_window_length);
- GI_sequence_left_position denotes the start and end positions of the predicted genomic islands.

### 5. Refining the boundaries of the predicted genomic islands using the MJSD-CG algorithm
- Find the max segmentation of a biological sequence using CG and MJD method
 GI_boundary(order_GI_predict,1)=CGMJD(test_seq1)+left_start;
 GI_boundary(order_GI_predict,2)=CGMJD(test_seq2)+right_start;
- GI_boundary denotes the refined start and end positions of the predicted genomic islands.



