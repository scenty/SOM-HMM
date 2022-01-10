clear 
%load MJK_res\CDOM\cdom_pattern.mat
%load MJK_res\som1020.mat %net
load D:\01-model_data\MJK\MJK_res\HMM\hmm_no2017
experiment_name='exp1'
load(['MJK_res\16n4\',experiment_name,'SOM.mat'])
val_label=1823:1989;
seq_com_cdom(t_model(:,2)>1989)=[];
type_cdom(t_model(:,2)>1989)=[];
t_model(t_model(:,2)>1989,:)=[];

[est_trans,est_emis]=hmmestimate(seq_com_model(train_label),type(train_label));
%training accuracy
likely_states=hmmviterbi(seq_com_model, est_trans, est_emis); 
p_train=sum(type(train_label)==likely_states(train_label))./length(train_label)
p_val =sum(type(val_label)==likely_states(val_label))./length(val_label)
%prediction accuracy
type_hmm=hmmviterbi(seq_com_cdom, est_trans, est_emis); 
p_predict=sum(type_cdom==type_hmm)./length(type_hmm)
%model accuracy
p_model=sum(type_cdom==type(t_model(:,2)))./length(type_cdom)

%save(['MJK_res\',experiment_name,'type_hmm.mat']) %type_hmm
save(['MJK_res\16n4\',experiment_name,'SOM.mat'],'type_hmm','p*','*label','-append')
