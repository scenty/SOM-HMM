clear 
load D:\01-model_data\MJK\MJK_res\HMM\cdom_shikong_no2017.mat sss_dailys_mean
experiment_name='exp9'
load(['MJK_res\16n4\',experiment_name,'SOM.mat'])

load MJK_res\mjk_lonlat.mat% mask
cdom_sss=sss_dailys_mean.*permute(repmat(mask,[1 1 size(sss_dailys_mean,1)]),[3 1 2]);


sk=5;cs=cdom_sss(:,1:sk:end,1:sk:end);
cs=reshape(cs,[size(cs,1),size(cs,2)*size(cs,3)]);
cs(:,isnan(sum(cs,1)))=[];
y_cdom=net(cs');
for iy=1:size(y_cdom,2)
    type_tmp(iy)=find(y_cdom(:,iy)==1);
    tmp=reshape([1:prod(som_neurons)],som_neurons)'; 
    tmp2=flipud(tmp);
    type_cdom(iy)=tmp2(find(type_tmp(iy)==tmp));
end

if 1
    type_tmp=type_cdom;
    type_tmp(type_cdom==4)=1;
    type_tmp(type_cdom==1)=3;
    type_tmp(type_cdom==3)=4;
    type_tmp(type_cdom==6)=5;
    type_tmp(type_cdom==5)=6;
    type_cdom=type_tmp;
end

save(['MJK_res\16n4\',experiment_name,'SOM.mat'],'type_cdom','-append')
