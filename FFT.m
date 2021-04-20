function out=FFT(data,fre,fs)
[~,fre_l]=size(fre);% 1*4
N=length(data);
fft_result=fft(data,N);  %每个点进行N次傅里叶，共N+1行，包括直流分量
modulus_value=abs(fft_result)*2/N; %求模，*2/N为固定公式
for i=1:fre_l    %%取各频率下的幅值和相角
        multiple=floor(fre(i)/fs*N)+1;%temp本质是倍频倍数，此处即为0 1 2 3
        %但0次分量在第1行，此处temp取1 2 3 4。floor为取整
        out(i,1)=modulus_value(multiple);
        out(i,2)=angle(fft_result(multiple));
    end
    out(1,1)=out(1,1)/2; %目前尚不清楚左右
end