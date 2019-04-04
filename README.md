# Topology_Inference_from_Consensus
Project 2: Network Topology Inference from Consensus Dynamics

#### Problem 1: Single consensus dynamics, known diffusion parameters
+ demo_p2_compare.m: compare InverseFilter, NearestCGL (NormF & Norm2), StructuredGLasso
  - compare the input is (a) white Gaussian (b) white uniform
  - compare the graph is (a) ER (36,0.1) (b) 6x6 grid (c) SBM (36,2,0.05,0.02) (d) smallworld (36,2,0.2) (e) scalefree (BA) (36,3,seed) the seed is set as a ER (5,0.2)
  
+ demo_p2_compare_exp.m  &  demo_p2_compare_betahop.m
  - compare algorithms under other filter functions: (a) exponential decay (b) power of L inverse

+ demo_p2_reason_consensus.m  &  demo_p2_reason_exp.m  &  demo_p2_reason_betahop.m
  - study the reason of performance differences under different filter functions

+ demo_p2_input_TimeSeries.m  &  demo_p2_input_Wishart.m
  - study the case when the input is not white

#### Problem 2: Single consensus dynamcis, constant and unknown diffusion parameters


#### Problem 3: Multiple consensus dynamcis
