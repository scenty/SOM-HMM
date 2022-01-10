clear;close all
ccode='v16';load(['MJK_res\',ccode,'SOM.mat'],'time')
%load MJK_res\som1020.mat %net

load D:\01-model_data\MJK\MJK_res\HMM\cdom_pattern_no2017.mat
load D:\01-model_data\MJK\MJK_res\HMM\hmm_no2017

%load D:\01-model_data\MJK\MJK_res\cdom.mat
%GSSS=32.5-28.7*cdom;


experiment_name='exp3' % 2*2, GvsH 2.61 MvH 1.92
experiment_name='exp1' % 3*2, GvsH 3.18 MvH 1.50 CV 1.41
%experiment_name='exp2' % 3*3, GvsH 2.56 MvH 1.54
%experiment_name='exp4' % 3*4, GvsH 2.7085 MvH 1.45
%experiment_name='exp5' % 4*4, GvsH 2.7080 MvH 1.37
%experiment_name='exp6' % 5*4, GvsH 2.69 MvH 1.28
%experiment_name='exp7' % 5*5, GvsH 2.70 MvH 1.28
%experiment_name='exp8' % 3*3,gridtop, GvsH 3.22 MvH 1.31 CV 1.20 
%experiment_name='exp9' % 3*3,  GvsH 3.21 MvH 1.3 CV 1.2
load(['MJK_res\16n4\',experiment_name,'SOM.mat'])
seq_com_cdom(t_model(:,2)>max(val_label))=[];
type_cdom(t_model(:,2)>max(val_label))=[];
t_model(t_model(:,2)>max(val_label),:)=[];
% estimate SSS with SOM
if ~exist('yhat','var')
for itype=1:prod(som_neurons)
   inx=find(type==itype);
   yhat(itype,:,:)=squeeze(nanmean(v1(inx,:,:)));
end
save(['MJK_res\16n4\',experiment_name,'SOM.mat'],'yhat','-append')
end
if ~exist('train_label','var')
    train_label=1:1500;val_label=1501:2002;
end
vhat(1:length(type),:,:)=yhat(type,:,:);
[rry,rsy]=drsGet2DR2(vhat,v1);
rrall=corrcoef(vhat(:),v1(:),'rows','complete').^2;
rsall=rmse(vhat(:),v1(:));


figure(1);clf
plot(time(t_model(:,2)),type(t_model(:,2)),'b.')
hold on
plot(t_model(:,1),type_cdom,'r.')
datetick
figure(2);clf
confusionchart(type_cdom,type(t_model(:,2)),'RowSummary','row-normalized');
confusionmat(type_cdom,type(t_model(:,2)))'/length(type_cdom)
title('GOCI vs Model')

figure(3);clf
confusionchart(type_cdom,type_hmm,'RowSummary','row-normalized');
confusionmat(type_cdom,type_hmm)'/length(type_hmm)
title('GOCI vs HMM')


figure(5);clf;set(gcf,'position',[1 1 560 420])
subplot(3,1,1)
plot(time,nanmean(vhat,[2 3]),'b')
hold on
plot(t_model(:,1),nanmean(cdom_sss,[2 3])-2,'r')
xlim([min(time) max(time)])
datetick('x','keeplimits')
tmp=vhat(t_model(:,2),:,:);
rsGvH=rmse(tmp(:),cdom_sss(:)-2);
title(['GOCI vs HMM RMSE=',num2str(rsGvH,'%.2f')])
subplot(3,1,2)
plot(time,nanmean(v1,[2 3]),'r')
hold on
plot(time,nanmean(vhat,[2 3]),'b')
xlim([min(time) max(time)])
datetick('x','keeplimits')
rsMvH=rmse(vhat(:),v1(:));
title(['Model vs HMM RMSE=',num2str(rsMvH,'%.2f')])
subplot(3,1,3)
plot(time,type)
xlim([min(time) max(time)])
datetick('x','keeplimits')
hold on
plot(time(type==2),type(type==2),'r.')

cff1=vhat(train_label,:,:);
cff2=v1(train_label,:,:);
rstrain=rmse(cff1(:),cff2(:));

cff1=vhat(val_label,:,:);
cff2=v1(val_label,:,:);
rsval=rmse(cff1(:),cff2(:));
figure(7);clf
pcolor(lon,lat,rsy)
shading interp
colorbar
caxis([0 4]);colormap jet
title(['Train RMSE=',num2str(rstrain,'%.2f'),' CV RMSE=',num2str(rsval,'%.2f')])

load D:\01-model_data\MJK\MJK_res\windspeed_direction1020.mat

for itype=1:length(unique(type))
umean(itype)=nanmean(U1_atime(type==itype));
vmean(itype)=nanmean(V1_atime(type==itype));
wmean(itype)=nanmean(w_atime(type==itype));
end

figure(8);clf
set(gcf,'color','w','position',[1111,428,500,450]);
boxplot(w_atime(1:length(type)),type,'Notch','on');
%ylabel('Wind Speed (m/s)','fontsize',16);
%xlabel('pattern','fontsize',16);
set(gca,'xticklabel',{'¢ñ','¢ò','¢ó','¢ô','¢õ','¢ö'},'yticklabel',[],'fontsize',16,'fontname','Î¢ÈíÑÅºÚ UI');
set(gca,'fontsize',18);
%myquiver(1:6,[1:6]*0+15,umean,vmean,'xref',4,'yref',12,'length',2,'scale',0.15,'arrowtype',.1,'linecolor',[1 0 1],'IsCenter')
grid on
ylim([0 20])

figure(9);clf
set(2,'color','w');
boxplot(V1_atime(1:length(type)),type,'Notch','on');
ylabel('V wind (m/s)','fontsize',16);
%xlabel('pattern','fontsize',16);
set(gca,'xticklabel',{'¢ñ','¢ò','¢ó','¢ô','¢õ','¢ö'},'fontsize',16);
set(gca,'fontsize',18);
set(9,'position',[1111,428,500,450]);

figure(10);clf;set(gcf,'position',[1 1 560 420])
subplot(2,1,1)
jj=41;ii=134;
plot(time,v1(:,jj,ii),'r')
hold on
plot(time,vhat(:,jj,ii),'b')
xlim([datenum(2020,1,1) datenum(2021,1,1)])
set(gca,'tickdir','out','linewidth',1,'fontsize',14,'xtick',datenum(2020,1:12,1))
datetick('x','m','keeplimits','keepticks')
rsMvH=rmse(v1(val_label,jj,ii),vhat(val_label,jj,ii));
title(['MJK Station CV-RMSE=',num2str(rsMvH,'%.2f'),' psu'])
ylabel('SSS (psu)')

subplot(2,1,2)
plot(time,nanmean(v1,[2 3]),'r')
hold on
plot(time,nanmean(vhat,[2 3]),'b')
xlim([datenum(2020,1,1) datenum(2021,1,1)])
set(gca,'tickdir','out','linewidth',1,'fontsize',14,'xtick',datenum(2020,1:12,1))
datetick('x','m','keeplimits','keepticks')
cff1=vhat(val_label,:,:);
cff2=v1(val_label,:,:);
rsMvH=rmse(cff1(:),cff2(:));
title(['Domain-averaged CV-RMSE=',num2str(rsMvH,'%.2f'),' psu'])
legend('Model SSS','HMM SSS')
ylabel('SSS (psu)')
xlabel('Month of 2020')


export_fig(figure(1),['mjk_res\type_Uwind_box'],'-jpg','-m2','-r2');
export_fig(figure(8),['mjk_res\type_Vwind_box'],'-tiff','-m3');
export_fig(figure(10),['mjk_res\PredictSSS'],'-tiff','-m3');

return

