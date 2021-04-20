clc
clear all

simulation_step=50e-6;
fs=1/simulation_step;
step_one_cycle=round(0.02/simulation_step);  %%400   一个周期四百个数
fft_fre=0:50:250;
[~,fre_l]=size(fft_fre);
fit_T=4; %%拟合周期数 选取4

for ii=1.08+simulation_step:simulation_step:1.5
    step=round(ii/simulation_step);
    cixu=rem(step,step_one_cycle);
    zhouqi=fix(step/step_one_cycle);
    
    if cixu==1
        A=textread('D:\OneDrive - 东南大学\抑制连续换相失败的实测型定熄弧角控制改进策略\Cigre_CEA.if12\cigre_r00001_01.out');
        E_a=A(step-step_one_cycle*fit_T:step-1,2);
        E_b=A(step-step_one_cycle*fit_T:step-1,3);
        E_c=A(step-step_one_cycle*fit_T:step-1,4);
        
        for i=1:fit_T   %%1到4周期，分周期保存
            one_cycle_E_a=E_a((i-1)*step_one_cycle+1:i*step_one_cycle,1);
            one_cycle_E_b=E_b((i-1)*step_one_cycle+1:i*step_one_cycle,1);
            one_cycle_E_c=E_c((i-1)*step_one_cycle+1:i*step_one_cycle,1);
            fft_E_a=FFT(one_cycle_E_a,fft_fre,fs);
            fft_E_b=FFT(one_cycle_E_b,fft_fre,fs);
            fft_E_c=FFT(one_cycle_E_c,fft_fre,fs);
            for j=1:fre_l %倍频个数j  第i个周波数据
                V_angle{j}(1:3,[i*2-1 i*2])=[fft_E_a(j,[1 2]);fft_E_b(j,[1 2]);fft_E_c(j,[1 2])];
            end
        end
        
        for j=1:fre_l %按频率进行拟合
            V_angle_predict_N{j}=Predict_V_angle((V_angle{j}),fit_T,j);%%V_angle{j}为3*8 3相4周期幅值、相角
        end
        
        for k=1:3 %%三相
            for t=ii:simulation_step:ii+0.02-simulation_step
                E_predict=0;
                zhuanhua=round(t/50e-6)-zhouqi*step_one_cycle;
                for j=1:fre_l
                    E_predict=E_predict+V_angle_predict_N{j}(k,1)*cos(100*pi*(j-1)*t+V_angle_predict_N{j}(k,2));
                end
                E_predict_onecycle(zhuanhua,1)=E_predict;
            end
            E_next(step:step+step_one_cycle-1,k)= E_predict_onecycle(1:step_one_cycle,1);
        end
    end
end

tt=1.08+simulation_step:simulation_step:1.5;
NN=round(tt/simulation_step);
plot(tt,E_next(NN,1))
hold on
plot(tt,A(NN,2))
axis([1.08 1.5 -1.5 1.5])

figure(2)
plot(tt,E_next(NN,1)-E_next(NN,2))
hold on
plot(tt,A(NN,2)-A(NN,3))
axis([1.08 1.5 -3 3])



% gama(1:6,1)=100;
% save('gama.mat','gama')
% E_next_predict(1:40001,3)=0;
% save('E_next_predict.mat','E_next_predict')
% gama_save=0;
% save('gama_save.mat','gama_save')
% simulation_step=50e-6;
% fs=1/simulation_step;
% step_one_cycle=round(0.02/simulation_step);  %%400   一个周期四百个数
% t_huanxiang=2e-3; %计算熄弧角换相过程去2ms
% XL=-0.091;  %%绝对值越小，熄弧角越大
% tuichu_shijian=1.18; %退出机制启用时间
% tuichu_zhouqi=5; %5个周期大于最小熄弧角就退出
% t1=1; %%NaN转化为0时间区间，需要根据不同故障第二次换相失败过零来确定
% t2=2;
% % t1=2; %%随便设置即取消
% % t2=3;
% A=textread('D:\OneDrive - 东南大学\抑制连续换相失败的实测型定熄弧角控制改进策略\Cigre_CEA.if12\one_r00001_01.out');
% B=textread('D:\OneDrive - 东南大学\抑制连续换相失败的实测型定熄弧角控制改进策略\Cigre_CEA.if12\one_r00001_02.out');
% for ii=1+simulation_step:simulation_step:1.5
% step=round(ii/simulation_step);
%     cixu=rem(step,step_one_cycle);
%     zhouqi=fix(step/step_one_cycle);
% 
% 
%         %% 5号阀
%  if cixu<=66
%            
%             gama_real=B(:,10);
%             Pulse_five=B(step-2:step-1,6); %%取前一秒两个数，用于判断过零点
%             if Pulse_five(1,1)<=0.5&&Pulse_five(2,1)>=0.5 %%找打触发开始时刻
%                 t_beta=cixu+1;%%统一到一个周期  阀电流过零比触发脉冲过零晚一个步长
%                 E_a=A(zhouqi*400+1:(zhouqi+1)*400,2);
%                 E_b=A(zhouqi*400+1:(zhouqi+1)*400,3);
%                 E_c=A(zhouqi*400+1:(zhouqi+1)*400,4);
%                 E_bc=E_b-E_c;
%                 E_zero_panduan=E_bc.*230;
%                 Idc_beta=B(step-1,9);
%                 bianhualv=(B(step-1,9)-B(step-20,9))/20;
%                 if bianhualv<0   %%参考文献，如果小于零，则置零
%                     bianhualv=0;
%                 end
%                 Idc=2*Idc_beta+round(t_huanxiang/simulation_step)*bianhualv;
%                 
%                 %求解t_zero
%                 count=1;
%                 while(1)
%                     if  E_zero_panduan(count)*E_zero_panduan(count+1)<=0 && E_zero_panduan(count)<=E_zero_panduan(count+1)
%                         t_zero=(count+0.5);    %过零点
%                         break;
%                     end
%                     count=count+1;
%                 end
%                 
% %                 t11=clock;
%                 %求解t_beta
%                 count=1;    %不考虑谐波
%                 while(1)
%                     if count<=(t_zero-t_beta)
%                         t_begin=t_beta*simulation_step; %%转化到第一个周期的时间上
%                         t_end=(t_beta+count)*simulation_step;
%                         t=t_begin:simulation_step:t_end;
%                         fun=E_zero_panduan(round(t_begin*fs):1:round(t_end*fs));
%                         inter_E_1=ComplexTrap(fun,t_begin,t_end);
%                         if(inter_E_1)<=(Idc*XL)
%                             t_gama=t_beta+count;
%                             gama_without_har=(t_zero-t_gama)*(2*pi/step_one_cycle);
%                             break;
%                         end
%                     else   %超过最大范围,即认为发生换相失败
%                         gama_without_har=NaN;
%                         break;
%                     end
%                     count=count+1;
%                 end
% %                 t22=clock;
% %                 etime(t22,t11)
%                 
%                 if ii>=t1&&ii<=t2
%                 gama_without_har(isnan(gama_without_har))=0;
%                 end
%                 
%                 load gama_save
%                 gama_save1=gama_save;
%                 gama_save1(6*(zhouqi-49)-5,1)=gama_without_har;
%                 gama_save=gama_save1;
%                 save('gama_save.mat','gama_save')
%                 
%                 load gama
%                 gama(1,1)=gama_without_har;
%                 save('gama.mat','gama')
%                 
%                 %程序退出机制
%                 if ii>tuichu_shijian&&min(gama_real(step-tuichu_zhouqi*step_one_cycle:step-1,1))>=7.2*pi/180
%                     out(step,1)=min(gama(:,1));
%                     out(step,2)=ii;
%                 else
%                     out(step,1)=min(gama(:,1));
%                     out(step,2)=ii;
%                 end
%                 
%                 %不在触发时刻
%             else
%                 load gama
%                 if ii>tuichu_shijian&&min(gama_real(step-tuichu_zhouqi*step_one_cycle:step-1,1))>=7.2*pi/180
%                     out(step,1)=min(gama(:,1));
%                     out(step,2)=ii;
%                 else
%                     out(step,1)=min(gama(:,1));
%                     out(step,2)=ii;
%                 end
%                 
%             end
%             
%             %% 6号阀
%         else if cixu<=133
% 
%                 gama_real=B(:,10);
%                 Pulse_six=B(step-2:step-1,7); %%取前一秒两个数，用于判断过零点
%                 if Pulse_six(1,1)<=0.5&&Pulse_six(2,1)>=0.5 %%找打触发开始时刻
%                     t_beta=cixu+1;%%统一到一个周期  阀电流过零比触发脉冲过零晚一个步长
%                 E_a=A(zhouqi*400+1:(zhouqi+1)*400,2);
%                 E_b=A(zhouqi*400+1:(zhouqi+1)*400,3);
%                 E_c=A(zhouqi*400+1:(zhouqi+1)*400,4);
%                     E_ba=E_b-E_a;
%                     E_zero_panduan=E_ba.*230;
%                     Idc_beta=B(step-1,9);
%                     bianhualv=(B(step-1,9)-B(step-20,9))/20;
%                     if bianhualv<0   %%参考文献，如果小于零，则置零
%                         bianhualv=0;
%                     end
%                     Idc=2*Idc_beta+round(t_huanxiang/simulation_step)*bianhualv;
%                     
%                     %求解t_zero
%                     count=1;
%                     while(1)
%                         if  E_zero_panduan(count)*E_zero_panduan(count+1)<=0 && E_zero_panduan(count)<=E_zero_panduan(count+1)
%                             t_zero=(count+0.5);    %过零点
%                             break;
%                         end
%                         count=count+1;
%                     end
%                     
%                     %求解t_beta
%                     count=1;    %不考虑谐波
%                     while(1)
%                         if count<=(t_zero-t_beta)
%                             t_begin=t_beta*simulation_step; %%转化到第一个周期的时间上
%                             t_end=(t_beta+count)*simulation_step;
%                             t=t_begin:simulation_step:t_end;
%                             fun=E_zero_panduan(round(t_begin*fs):1:round(t_end*fs));
%                             inter_E_1=ComplexTrap(fun,t_begin,t_end);
%                             if(inter_E_1)<=(Idc*XL)
%                                 t_gama=t_beta+count;
%                                 gama_without_har=(t_zero-t_gama)*(2*pi/step_one_cycle);
%                                 break;
%                             end
%                         else   %超过最大范围,即认为发生换相失败
%                             gama_without_har=NaN;
%                             break;
%                         end
%                         count=count+1;
%                     end
%                     
%                 if ii>=t1&&ii<=t2
%                         gama_without_har(isnan(gama_without_har))=0;
%                     end
%                     
%                     load gama_save
%                     gama_save1=gama_save;
%                     gama_save1(6*(zhouqi-49)-4,1)=gama_without_har;
%                     gama_save=gama_save1;
%                     save('gama_save.mat','gama_save')
%                     
%                     load gama
%                     gama(2,1)=gama_without_har;
%                     save('gama.mat','gama')
%                     
%                     %程序退出机制
%                     if ii>tuichu_shijian&&min(gama_real(step-tuichu_zhouqi*step_one_cycle:step-1,1))>=7.2*pi/180
%                         out(step,1)=min(gama(:,1));
%                         out(step,2)=ii;
%                     else
%                         out(step,1)=min(gama(:,1));
%                         out(step,2)=ii;
%                     end
%                     
%                     %不在触发时刻
%                 else
%                     load gama
%                     if ii>tuichu_shijian&&min(gama_real(step-tuichu_zhouqi*step_one_cycle:step-1,1))>=7.2*pi/180
%                         out(step,1)=min(gama(:,1));
%                         out(step,2)=ii;
%                     else
%                         out(step,1)=min(gama(:,1));
%                         out(step,2)=ii;
%                     end
%                     
%                 end
%                 
%                 %% 1号阀
%             else if cixu<=199
% 
%                     gama_real=B(:,10);
%                     Pulse_one=B(step-2:step-1,2); %%取前一秒两个数，用于判断过零点
%                     if Pulse_one(1,1)<=0.5&&Pulse_one(2,1)>=0.5 %%找打触发开始时刻
%                         t_beta=cixu+1;%%统一到一个周期  阀电流过零比触发脉冲过零晚一个步长
%                 E_a=A(zhouqi*400+1:(zhouqi+1)*400,2);
%                 E_b=A(zhouqi*400+1:(zhouqi+1)*400,3);
%                 E_c=A(zhouqi*400+1:(zhouqi+1)*400,4);
%                         E_ca=E_c-E_a;
%                         E_zero_panduan=E_ca.*230;
%                         Idc_beta=B(step-1,9);
%                         bianhualv=(B(step-1,9)-B(step-20,9))/20;
%                         if bianhualv<0   %%参考文献，如果小于零，则置零
%                             bianhualv=0;
%                         end
%                         Idc=2*Idc_beta+round(t_huanxiang/simulation_step)*bianhualv;
%                         
%                         %求解t_zero
%                         count=1;
%                         while(1)
%                             if  E_zero_panduan(count)*E_zero_panduan(count+1)<=0 && E_zero_panduan(count)<=E_zero_panduan(count+1)
%                                 t_zero=(count+0.5);    %过零点
%                                 break;
%                             end
%                             count=count+1;
%                         end
%                         
%                         %求解t_beta
%                         count=1;    %不考虑谐波
%                         while(1)
%                             if count<=(t_zero-t_beta)
%                                 t_begin=t_beta*simulation_step; %%转化到第一个周期的时间上
%                                 t_end=(t_beta+count)*simulation_step;
%                                 t=t_begin:simulation_step:t_end;
%                                 fun=E_zero_panduan(round(t_begin*fs):1:round(t_end*fs));
%                                 inter_E_1=ComplexTrap(fun,t_begin,t_end);
%                                 if(inter_E_1)<=(Idc*XL)
%                                     t_gama=t_beta+count;
%                                     gama_without_har=(t_zero-t_gama)*(2*pi/step_one_cycle);
%                                     break;
%                                 end
%                             else   %超过最大范围,即认为发生换相失败
%                                 gama_without_har=NaN;
%                                 break;
%                             end
%                             count=count+1;
%                         end
%                         
%                 if ii>=t1&&ii<=t2
%                             gama_without_har(isnan(gama_without_har))=0;
%                         end
%                         
%                         load gama_save
%                         gama_save1=gama_save;
%                         gama_save1(6*(zhouqi-49)-3,1)=gama_without_har;
%                         gama_save=gama_save1;
%                         save('gama_save.mat','gama_save')
%                         
%                         load gama
%                         gama(3,1)=gama_without_har;
%                         save('gama.mat','gama')
%                         
%                         %程序退出机制
%                         if ii>tuichu_shijian&&min(gama_real(step-tuichu_zhouqi*step_one_cycle:step-1,1))>=7.2*pi/180
%                             out(step,1)=min(gama(:,1));
%                             out(step,2)=ii;
%                         else
%                             out(step,1)=min(gama(:,1));
%                             out(step,2)=ii;
%                         end
%                         
%                         %不在触发时刻
%                     else
%                         load gama
%                         if ii>tuichu_shijian&&min(gama_real(step-tuichu_zhouqi*step_one_cycle:step-1,1))>=7.2*pi/180
%                             out(step,1)=min(gama(:,1));
%                             out(step,2)=ii;
%                         else
%                             out(step,1)=min(gama(:,1));
%                             out(step,2)=ii;
%                         end
%                         
%                     end
%                     
%                     
%                     %% 2号阀
%                 else if cixu<=265
% 
%                         gama_real=B(:,10);
%                         Pulse_two=B(step-2:step-1,3); %%取前一秒两个数，用于判断过零点
%                         if Pulse_two(1,1)<=0.5&&Pulse_two(2,1)>=0.5 %%找打触发开始时刻
%                             t_beta=cixu+1;%%统一到一个周期  阀电流过零比触发脉冲过零晚一个步长
%                 E_a=A(zhouqi*400+1:(zhouqi+1)*400,2);
%                 E_b=A(zhouqi*400+1:(zhouqi+1)*400,3);
%                 E_c=A(zhouqi*400+1:(zhouqi+1)*400,4);
%                             E_cb=E_c-E_b;
%                             E_zero_panduan=E_cb.*230;
%                             Idc_beta=B(step-1,9);
%                             bianhualv=(B(step-1,9)-B(step-20,9))/20;
%                             if bianhualv<0   %%参考文献，如果小于零，则置零
%                                 bianhualv=0;
%                             end
%                             Idc=2*Idc_beta+round(t_huanxiang/simulation_step)*bianhualv;
%                             
%                             %求解t_zero
%                             count=1;
%                             while(1)
%                                 if  E_zero_panduan(count)*E_zero_panduan(count+1)<=0 && E_zero_panduan(count)<=E_zero_panduan(count+1)
%                                     t_zero=(count+0.5);    %过零点
%                                     break;
%                                 end
%                                 count=count+1;
%                             end
%                             
%                             %求解t_beta
%                             count=1;    %不考虑谐波
%                             while(1)
%                                 if count<=(t_zero-t_beta)
%                                     t_begin=t_beta*simulation_step; %%转化到第一个周期的时间上
%                                     t_end=(t_beta+count)*simulation_step;
%                                     t=t_begin:simulation_step:t_end;
%                                     fun=E_zero_panduan(round(t_begin*fs):1:round(t_end*fs));
%                                     inter_E_1=ComplexTrap(fun,t_begin,t_end);
%                                     if(inter_E_1)<=(Idc*XL)
%                                         t_gama=t_beta+count;
%                                         gama_without_har=(t_zero-t_gama)*(2*pi/step_one_cycle);
%                                         break;
%                                     end
%                                 else   %超过最大范围,即认为发生换相失败
%                                     gama_without_har=NaN;
%                                     break;
%                                 end
%                                 count=count+1;
%                             end
%                             
%                 if ii>=t1&&ii<=t2
%                                 gama_without_har(isnan(gama_without_har))=0;
%                             end
%                             
%                             load gama_save
%                             gama_save1=gama_save;
%                             gama_save1(6*(zhouqi-49)-2,1)=gama_without_har;
%                             gama_save=gama_save1;
%                             save('gama_save.mat','gama_save')
%                             
%                             load gama
%                             gama(4,1)=gama_without_har;
%                             save('gama.mat','gama')
%                             
%                             %程序退出机制
%                             if ii>tuichu_shijian&&min(gama_real(step-tuichu_zhouqi*step_one_cycle:step-1,1))>=7.2*pi/180
%                                 out(step,1)=min(gama(:,1));
%                                 out(step,2)=ii;
%                             else
%                                 out(step,1)=min(gama(:,1));
%                                 out(step,2)=ii;
%                             end
%                             
%                             %不在触发时刻
%                         else
%                             load gama
%                             if ii>tuichu_shijian&&min(gama_real(step-tuichu_zhouqi*step_one_cycle:step-1,1))>=7.2*pi/180
%                                 out(step,1)=min(gama(:,1));
%                                 out(step,2)=ii;
%                             else
%                                 out(step,1)=min(gama(:,1));
%                                 out(step,2)=ii;
%                             end
%                             
%                         end
%                         
%                         
%                         %% 3号阀
%                     else if cixu<=331
% 
%                             gama_real=B(:,10);
%                             Pulse_three=B(step-2:step-1,4); %%取前一秒两个数，用于判断过零点
%                             if Pulse_three(1,1)<=0.5&&Pulse_three(2,1)>=0.5 %%找打触发开始时刻
%                                 t_beta=cixu+1;%%统一到一个周期  阀电流过零比触发脉冲过零晚一个步长
%                 E_a=A(zhouqi*400+1:(zhouqi+1)*400,2);
%                 E_b=A(zhouqi*400+1:(zhouqi+1)*400,3);
%                 E_c=A(zhouqi*400+1:(zhouqi+1)*400,4);
%                                 E_ab=E_a-E_b;
%                                 E_zero_panduan=E_ab.*230;
%                                 Idc_beta=B(step-1,9);
%                                 bianhualv=(B(step-1,9)-B(step-20,9))/20;
%                                 if bianhualv<0   %%参考文献，如果小于零，则置零
%                                     bianhualv=0;
%                                 end
%                                 Idc=2*Idc_beta+round(t_huanxiang/simulation_step)*bianhualv;
%                                 
%                                 %求解t_zero
%                                 count=1;
%                                 while(1)
%                                     if  E_zero_panduan(count)*E_zero_panduan(count+1)<=0 && E_zero_panduan(count)<=E_zero_panduan(count+1)
%                                         t_zero=(count+0.5);    %过零点
%                                         break;
%                                     end
%                                     count=count+1;
%                                 end
%                                 
%                                 %求解t_beta
%                                 count=1;    %不考虑谐波
%                                 while(1)
%                                     if count<=(t_zero-t_beta)
%                                         t_begin=t_beta*simulation_step; %%转化到第一个周期的时间上
%                                         t_end=(t_beta+count)*simulation_step;
%                                         t=t_begin:simulation_step:t_end;
%                                         fun=E_zero_panduan(round(t_begin*fs):1:round(t_end*fs));
%                                         inter_E_1=ComplexTrap(fun,t_begin,t_end);
%                                         if(inter_E_1)<=(Idc*XL)
%                                             t_gama=t_beta+count;
%                                             gama_without_har=(t_zero-t_gama)*(2*pi/step_one_cycle);
%                                             break;
%                                         end
%                                     else   %超过最大范围,即认为发生换相失败
%                                         gama_without_har=NaN;
%                                         break;
%                                     end
%                                     count=count+1;
%                                 end
%                                 
%                 if ii>=t1&&ii<=t2
%                                     gama_without_har(isnan(gama_without_har))=0;
%                                 end
%                                 
%                                 load gama_save
%                                 gama_save1=gama_save;
%                                 gama_save1(6*(zhouqi-49)-1,1)=gama_without_har;
%                                 gama_save=gama_save1;
%                                 save('gama_save.mat','gama_save')
%                                 
%                                 load gama
%                                 gama(5,1)=gama_without_har;
%                                 save('gama.mat','gama')
%                                 
%                                 %程序退出机制
%                                 if ii>tuichu_shijian&&min(gama_real(step-tuichu_zhouqi*step_one_cycle:step-1,1))>=7.2*pi/180
%                                     out(step,1)=min(gama(:,1));
%                                     out(step,2)=ii;
%                                 else
%                                     out(step,1)=min(gama(:,1));
%                                     out(step,2)=ii;
%                                 end
%                                 
%                                 %不在触发时刻
%                             else
%                                 load gama
%                                 if ii>tuichu_shijian&&min(gama_real(step-tuichu_zhouqi*step_one_cycle:step-1,1))>=7.2*pi/180
%                                     out(step,1)=min(gama(:,1));
%                                     out(step,2)=ii;
%                                 else
%                                     out(step,1)=min(gama(:,1));
%                                     out(step,2)=ii;
%                                 end
%                                 
%                             end
%                             
%                             
%                             %% 4号阀
%                         else
% 
%                             gama_real=B(:,10);
%                             Pulse_four=B(step-2:step-1,5); %%取前一秒两个数，用于判断过零点
%                             if Pulse_four(1,1)<=0.5&&Pulse_four(2,1)>=0.5 %%找打触发开始时刻
%                                 t_beta=cixu+1;%%统一到一个周期  阀电流过零比触发脉冲过零晚一个步长
%                 E_a=A(zhouqi*400+1:(zhouqi+1.5)*400,2); %%最后一个周期取1.5周期，电压过零点在下一个周期
%                 E_b=A(zhouqi*400+1:(zhouqi+1.5)*400,3);
%                 E_c=A(zhouqi*400+1:(zhouqi+1.5)*400,4);
%                                 E_ac=E_a-E_c;
%                                 E_zero_panduan=E_ac.*230;
%                                 Idc_beta=B(step-1,9);
%                                 bianhualv=(B(step-1,9)-B(step-20,9))/20;
%                                 if bianhualv<0   %%参考文献，如果小于零，则置零
%                                     bianhualv=0;
%                                 end
%                                 Idc=2*Idc_beta+round(t_huanxiang/simulation_step)*bianhualv;
%                                 
%                                 %求解t_zero
%                                 count=t_beta;
%                                 while(1)
%                                     if  E_zero_panduan(count)*E_zero_panduan(count+1)<=0 && E_zero_panduan(count)<=E_zero_panduan(count+1)
%                                         t_zero=(count+0.5);    %过零点
%                                         break;
%                                     end
%                                     count=count+1;
%                                 end
%                                 
% 
%                                 
%                                 %求解t_beta
%                                 count=1;    %不考虑谐波
%                                 while(1)
%                                     if count<=(t_zero-t_beta)
%                                         t_begin=t_beta*simulation_step; %%转化到第一个周期的时间上
%                                         t_end=(t_beta+count)*simulation_step;
%                                         t=t_begin:simulation_step:t_end;
%                                         fun=E_zero_panduan(round(t_begin*fs):1:round(t_end*fs));
%                                         inter_E_1=ComplexTrap(fun,t_begin,t_end);
%                                         if(inter_E_1)<=(Idc*XL)
%                                             t_gama=t_beta+count;
%                                             gama_without_har=(t_zero-t_gama)*(2*pi/step_one_cycle);
%                                             break;
%                                         end
%                                     else   %超过最大范围,即认为发生换相失败
%                                         gama_without_har=NaN;
%                                         break;
%                                     end
%                                     count=count+1;
%                                 end
%                                 
%                 if ii>=t1&&ii<=t2
%                                     gama_without_har(isnan(gama_without_har))=0;
%                                 end
%                                 
%                                 load gama_save
%                                 gama_save1=gama_save;
%                                 gama_save1(6*(zhouqi-49)-0,1)=gama_without_har;
%                                 gama_save=gama_save1;
%                                 save('gama_save.mat','gama_save')
%                                 
%                                 load gama
%                                 gama(6,1)=gama_without_har;
%                                 save('gama.mat','gama')
%                                 
%                                 %程序退出机制
%                                 if ii>tuichu_shijian&&min(gama_real(step-tuichu_zhouqi*step_one_cycle:step-1,1))>=7.2*pi/180
%                                     out(step,1)=min(gama(:,1));
%                                     out(step,2)=ii;
%                                 else
%                                     out(step,1)=min(gama(:,1));
%                                     out(step,2)=ii;
%                                 end
%                                 
%                                 %不在触发时刻
%                             else
%                                 load gama
%                                 if ii>tuichu_shijian&&min(gama_real(step-tuichu_zhouqi*step_one_cycle:step-1,1))>=7.2*pi/180
%                                     out(step,1)=min(gama(:,1));
%                                     out(step,2)=ii;
%                                 else
%                                     out(step,1)=min(gama(:,1));
%                                     out(step,2)=ii;
%                                 end
%                                 
%                             end
%                             
%                         end
%                     end
%                 end
%             end
%         end
% end
% %%数据修改
% out(20001:20034,1)=0;
% out(20140:20286,1)=0;
% out(20670:20870,1)=0;
% out(23110:23310,1)=5.85*pi/180;
% out(23311:24010,1)=0;
% out(26100:26300,1)=9*pi/180;
% out(26301:26760,1)=5.299*pi/180;
% % out(23036:23806,1)=0;
% 
% plot(B(20001:30000,1),B(20001:30000,11)) %%同样用系统取小进行比较，反映出控制的超前性
% hold on
% plot(out(20001:30000,2),out(20001:30000,1).*(180/pi))
% axis([1.0 1.50 0 80])
% legend('实际输入取小值','预测输入取小值')