%close all;clear all;
swat=load('D:\01-model_data\MJK\runoff\SWAT_daily.mat');
%
experiment_name='exp8'
plot_tracer=1;
plot_map=0; %overlap google map or not?
plot_zeta=0;
plot_currents=0;
from_nc=0;
export_map=1; % experot or not?
resample = 1;
train_label=1:1822;
fdir='D:\01-model_data\MJK\res\';
%files=dir([fdir,'qck.mjk.v16*.nc']);
files=dir([fdir,'qck.mjk.v16*.nc']);
lon=nc_varget('D:\01-model_data\MJK\res\qck.mjk.v16.nc','lon_rho');
lat=nc_varget('D:\01-model_data\MJK\res\qck.mjk.v16.nc','lat_rho');
ccode=chooseCase('D:\01-model_data\MJK\res\qck.mjk.v16.nc');
if from_nc
    [time,atime,v0,z0]=deal([]);
    for iff=1:length(files)
        fname=[fdir,files(iff).name];
        aname=fname;ik=strfind(aname,'qck');
        aname(ik:ik+2)='avg';
        time=[time; nc_varget(aname,'ocean_time')/86400+datenum(1900,1,1)];
        %atime=[atime; nc_varget(aname,'ocean_time')/86400+datenum(1900,1,1)];
        %v0=[v0; nc_varget(fname,'salt_sur')];
        %z0=[z0; nc_varget(fname,'zeta')];
        nc=netcdf(aname);    tmp=squeeze(nc{'salt'}(:,10,:,:));   close(nc)
        v0=[v0; tmp];
        z0=[z0; nc_varget(aname,'zeta')];
    end
    v0(v0>1e20)=NaN;z0(z0>1e20)=NaN;
    
    %figure(1)
    %pcolor(lon,lat,squeeze(nanmean(v0)))
    %shading interp;colorbar
    if resample==1
        %reduce v1 [N*M] to [2*N/2*M]
        dt = 3;
        [N,M]=size(v0);
        sm=mod(N,dt);
        time=time(1:end-sm); v0=v0(1:end-sm,:,:); z0=z0(1:end-sm,:,:);
        time=nanmean(reshape(time,[dt length(time)/dt]));
        v1=squeeze(nanmean(reshape(v0,[dt length(v0)/dt size(v0,2) size(v0,3)])));
        z1=squeeze(nanmean(reshape(z0,[dt length(z0)/dt size(z0,2) size(z0,3)])));
        %%%%%%caution!
    elseif resample==2
        dt=2;
            for ii=1:length(atime)%（特别重要：时间的匹配方式）
               inx = find(time>= atime(ii)-dt/2 & time<= atime(ii)+dt/2 );
               nc=netcdf(fname);
               tmp=v0(inx,:,:);
               %disp(inx)
               if length(inx)~=3; error(''); end
               v1(ii,:,:)=nanmean(tmp);
            end
    else
        v1=v0;clear v0
    end
else
    load D:\01-model_data\MJK\MJK_res\ModelSSS2010-2020.16n4.mat
    s_sur(s_sur>1e10)=NaN;
    time=datenum(2010,1,1:2:2*length(s_sur));
end
v1=s_sur;
sk=5;
v2=v1(:,1:sk:end,1:sk:end);%reducing the dimension
xx=lon(1:sk:end,1:sk:end);yy=lat(1:sk:end,1:sk:end);

%v2=v2(:,find(xx<=119.75));
v2=reshape(v2,[size(v2,1) size(v2,2)*size(v2,3)]);
mask= (~isnan(squeeze((sum(v1(:,1:sk:end,1:sk:end),1)))) & xx<=119.75);
v2(:,isnan(sum(v2,1)))=[];
% figure(2);clf
% pcolor(lon1,lat1,squeeze(nanmean(v1(:,1:sk:end,1:sk:end))))
% shading interp
%interpolation of flow
flow=interp1(swat.dt,swat.flow_daily,time);
%
som_neurons=[3 3];
net = selforgmap(som_neurons,100,3,'gridtop');
[net,tr] = train(net,v2(train_label,:)');
%
y=net(v2');
%type=mod(find(y==1),16);

%%%%%%%%%%
tmp=reshape([1:prod(som_neurons)],som_neurons)';
tmp2=flipud(tmp); %important
for iy=1:size(y,2)
type_tmp(iy)=find(y(:,iy)==1);
%or type_tmp = vec2ind(net(v1')); by LWF
%[13 14 15 16]
%[ 9 10 11 12]
%[ 5  6  7  8]
%[ 1  2  3  4]
%inx=find(type==tmp2(itype)); %right order
type(iy)=tmp2(find(type_tmp(iy)==tmp));
end
%%%%%%%%%%
if 1
    type_tmp=type;
    type_tmp(type==4)=1;
    type_tmp(type==1)=3;
    type_tmp(type==3)=4;
    type_tmp(type==6)=5;
    type_tmp(type==5)=6;
    type=type_tmp;
end

%load([ccode,'SOM.mat'])

figure(4);clf;set(gcf,'position',[1 1 928 800])
s=tight_subplot(som_neurons(2),som_neurons(1),0.01);
for itype=1:prod(som_neurons)
   
   inx=find(type==itype);
   subplot(s(itype))
   yhat(itype,:,:)=squeeze(nanmean(v1(inx,:,:)));
   pcolor(lon,lat,squeeze(nanmean(v1(inx,:,:))));
   text(119.47,26.3,{['flow=',num2str(nanmean(flow(inx))/1e4,'%.2f'),'\pm',num2str(nanstd(flow(inx))/1e4,'%.2f')],[num2str(length(inx)./length(type)*100,'%3.1f'),'%']},'fontsize',12)
   hold on
   %contour(lon,lat,squeeze(nanmean(v0(inx,:,:))),[20 20],'k-');
   if plot_zeta
   mz=nanmean(squeeze(nanmean(z1(inx,:,:))),[1 2]);
   contour(lon,lat,squeeze(nanmean(z1(inx,:,:)))-mz,[0:.01:.05],'k-');
   contour(lon,lat,squeeze(nanmean(z1(inx,:,:)))-mz,[-.05:.01:0],'k--');
   contour(lon,lat,mz-squeeze(nanmean(z1(inx,:,:))),[0 0],'k-','linewidth',2);
   end
   %clabel(cl,ch,[.9 .9],'labelspacing',100)
   shading interp
   caxis([14 34])
   cptcmap('haline','ncol',20)
    if itype<=som_neurons(1)*(som_neurons(1)-1)
        set(gca,'xtick',[])
    end
    if mod(itype-1,som_neurons(1))~=0
        set(gca,'ytick',[])
    end
    set(gca,'fontsize',12)
end
save(['mjk_res\16n4\',experiment_name,'SOM.mat'],'time','type','mask','lon','lat','v1','som_neurons','net','yhat','train_label')

vhat(1:length(type),:,:)=yhat(type,:,:);

figure;set(gcf,'position',[1600 154 560 420])
plotsomnd(net,v2')
figure
plotsomhits(net,v2')

%mjk_compositeSec
return
lag=7;
clear flow_feature
N2=length(flow)-lag
for ii=1:N2
    flow_feature(ii,:)=flow(ii:ii+lag-1);
end
net2 = selforgmap([2 2]);
[net2,tr2] = train(net2,flow_feature');
%
%y=net2(flow_feature');
type_flow = vec2ind(net2(flow_feature'));

figure(6);clf;set(gcf,'position',[1600 600 560 420])
subplot(2,1,1)
plot(time(1:N2),type_flow)
xlim([min(time) max(time)])
datetick('x','keeplimits')
subplot(2,1,2)
plot(time(1:N2),flow(1:N2))
xlim([min(time) max(time)])
datetick('x','keeplimits')
figure
plotsompos(net2)
