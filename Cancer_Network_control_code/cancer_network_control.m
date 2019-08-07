function [ MMS,MDS,DFVS,NCUA ] = cancer_network_control( expression_tumor_fileName,expression_normal_fileName,Network_method_index )
%we output the sample-specific driver profiles by using different control
%methods
%   Input:
%         expression_fileName including expression_tumor_fileName and expression_normal_fileName)   
%         index:denotes we use which network construction method
%**************the network information****
%if index=1,we use SSN
%if Network_index=2,we use paired SSN; 
%if Network_index=3,we use LIONESS
%  Output:
%         The sample-specific driver profiles of MMS,MDS,DFVS,NCUA;
%         The column is the samples and the rows is the genes. The value “1” denoted that the gene is driver genes; 
%************************part1:LOAD sample data and network data************************
%********************obtain the paired expression data******************

[tumor,~,name_tumor]=importdata(expression_tumor_fileName);
gene_list=tumor.textdata(2:end,1);tumor_data=tumor.data;
[normal,~,name_normal]=importdata(expression_normal_fileName);
Sample_name_normal=normal.textdata(1,2:end);normal_data=normal.data;
data=tumor_data;ref_data=normal_data;


index=Network_method_index;
load('GIN_network_information.mat')
%load('sPPI_network_information.mat')
%[~,PPI,~]=xlsread('network_FIsInGene_041709.xlsx');
[x1,y1]=ismember(edge0(:,1),gene_list);
[x2,y2]=ismember(edge0(:,2),gene_list);

y=y1.*y2;
z=[y1 y2];
z(find(y==0),:)=[];
N1=length(gene_list);
[N2,~]=size(z);

Net=zeros(N1,N1);

for i=1:N2
    
         Net(z(i,2),z(i,1))=1;  %undirected gene-gene interaction network
         %Net(z(i,1),z(i,2))=1;    
end



if index==1
    
  
%***************SSN*********************** 
 
for i=1:size(data,2)
%  for i=1:2
      
    tic
    i
    sample=data(:,i);
    [ index_R,p ] =SSN(sample,ref_data);
    p(isnan(p))=1;
    p(p>=0.05)=0;
    p(p~=0)=1;
    
    cand=p.*Net;
    
    C=cand;
    [x,y]=find(C~=0);
    Dz1=[y x];N=length(C);

    
[ MMS_x1,MMS_nd1 ] = control( Dz1,N );
[ MDS_x1,MDS_nd1 ] = Opti_MDS( Dz1,N  );
[ NCD_x1,NCD_nd1 ] = Opti_weight_ncd( Dz1,N  );
[ NCU_x1,NCU_nd1 ] = Opti_weight_nc( Dz1,N  );

MMS(:,i)=MMS_x1;
MDS(:,i)=MDS_x1;
DFVS(:,i)=NCD_x1;
NCUA(:,i)=NCU_x1;


    clear index_R p cand

 toc   
end




end



%ndm = csndm(a,0.01,0.1,0);

if index==2
    
 %***************paired SSN*********************** 
 
  
for i=1:size(data,2)
%  for i=1:2
    
    tic
    i
   
    sample_tumor=data(:,i);
    [R0,P]=SSN(sample_tumor,ref_data);
    
   P(isnan(P))=1;
   P(P>=0.05)=0;
   P(P~=0)=1;
    %construct the normal SSN
    clear  sample_tumor 
    sample_normal=ref_data(:,i);
    [R1,P1]=SSN(sample_normal,ref_data);
    clear  sample_normal 
    
   P1(isnan(P1))=1;
   P1(P1>=0.05)=0;
   P1(P1~=0)=1;
    C=abs(P-P1).*Net; 
    
    
    [x,y]=find(C~=0);
    Dz1=[y x];N=length(C);

    
[ MMS_x1,MMS_nd1 ] = control( Dz1,N );
[ MDS_x1,MDS_nd1 ] = Opti_MDS( Dz1,N  );
[ NCD_x1,NCD_nd1 ] = Opti_weight_ncd( Dz1,N  );
[ NCU_x1,NCU_nd1 ] = Opti_weight_nc( Dz1,N  );

MMS(:,i)=MMS_x1;
MDS(:,i)=MDS_x1;
DFVS(:,i)=NCD_x1;
NCUA(:,i)=NCU_x1;
    
 toc   
 
    
end



end




if index==3
    
%***************LIONESS***********************   
  
for i=1:size(data,2)
%  for i=1:2
     tic   
   i
    p=lioness_method(data,i);
  
    C=p.*Net; 

    [x,y]=find(C~=0);
    Dz1=[y x];N=length(C);

   
[ MMS_x1,MMS_nd1 ] = control( Dz1,N );
[ MDS_x1,MDS_nd1 ] = Opti_MDS( Dz1,N  );
[ NCD_x1,NCD_nd1 ] = Opti_weight_ncd( Dz1,N  );
[ NCU_x1,NCU_nd1 ] = Opti_weight_nc( Dz1,N  );

MMS(:,i)=MMS_x1;
MDS(:,i)=MDS_x1;
DFVS(:,i)=NCD_x1;
NCUA(:,i)=NCU_x1;
    
clear C
    
   toc
    
   end
   
   

   
end




end

