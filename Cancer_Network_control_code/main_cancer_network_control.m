clc
clear
%   $Id: main_Benchmark_control.m Created at 2019-05-08 18:15:16 $
%   by Weifeng Guo, Northwestern Polytechtical University, China
%   Copyright (c) 2014-2019 by Key Laboratory of Information Fusion Technology of Ministry of Education in Northwestern Polytechnical University,
%   and key Laboratory of Systems Biology in Shanghai Institutes for Biological Science; 
%   If any problem,pleasse contact shaonianweifeng@126.com for help.

%Remainder: Please install gurobi before running our code
%Remainder: Please install gurobi before running our code
%Remainder: Please install gurobi before running our code

%**************Part 1:Input the information of samples and network information****
%**************sample information**************
%Example:a demo containing data of 5 patients
expression_tumor_fileName = 'Example_tumor.txt';
expression_normal_fileName = 'Example_normal.txt';


%**************the network construction information****
%if Network_index=1,we use SSN; if Network_index=2,we use paired SSN
%if Network_index=3,we use LIONESS
Network_method_index=1;
%Network_method_index=2;
%Network_method_index=3;

%%**************Part 2:Network control methods output the predicted combinational drugs****

[ MMS,MDS,DFVS,NCUA ] = cancer_network_control( expression_tumor_fileName,expression_normal_fileName,Network_method_index );

%%**************Part 3:save the result****

save Cancer_network_control_results
