%%From the response UA of a VPL neuron to a stimulus the code will
%%calculate the FR and the dynamic non linear prediction (as well as VAF of
%%the fit and ambiguity value).

warning off;
%%Load the response UA of a VPL neuron to a stimulus  
load stimulus.mat
load UA.mat

T1=0.06; T2=0.0006; Tc=5.7; k=1.7207;
Pars_transfer=[T1 Tc T2 k]; Pars_amplitude=[2 0.037];

nfft=2048;
noverlap=1024;
dflag='none';
n_tapers=8;
[E,V]=dpss(2048,n_tapers/2);

t=1+1/1000*(1:length(stimulus));
t=t(:);
time=t;
H=hilbert(stimulus-mean(stimulus));
amplitude=abs(H);
phase_component=phase(H);
frequency=phase_component./2./pi./t;

[pxx,fyy]=pwelch(stimulus,bartlett(2048),1024,2048,1000);
transfer_function=transfer(Pars_transfer,fyy);
gain=abs(transfer_function);
phase_transfer=phase(transfer_function);
transfer_function_neg=conj(transfer_function(length(transfer_function):-1:1));
filter=real(ifft([transfer_function' transfer_function_neg']'));
filter=filter([length(filter)/2:length(filter) 1:length(filter)/2-1]);
[b,a]=butter(1,15/500);
filter_filt=filtfilt(b,a,filter);
transformed_stimulus_linear=conv(stimulus,filter_filt,'same');
Gain_linear=gain;
phase_linear=phase_transfer*180/pi;
figure;subplot(2,1,1);plot(fyy,Gain_linear);subplot(2,1,2);plot(fyy,phase_linear);
title('linear')

transformed_stimulus=abs(transfer(Pars_transfer,frequency)).*amplitude.*(0+1*fitfunction(Pars_amplitude,amplitude)).*cos(2*pi.*frequency.*t+atan2(imag(transfer(Pars_transfer,frequency)),real(transfer(Pars_transfer,frequency))));
transformed_stimulus_new=sigmoid([-25 50 0 10],transformed_stimulus);
transformed_stimulus=transformed_stimulus_new;

[firing_rate, kaiser_filt] = kaiser_filter_lowpass(20,UA,1000);
predicted_firing_rate_linear=transformed_stimulus;%transformed_stimulus_linear
bins=linspace(min(predicted_firing_rate_linear),max(predicted_firing_rate_linear),100);
for I=1:length(bins)-1
    mean_rate(I)=median(firing_rate(find(predicted_firing_rate_linear<=bins(I+1) & predicted_firing_rate_linear>bins(I))));
    ste_rate(I)=std(firing_rate(find(predicted_firing_rate_linear<=bins(I+1) & predicted_firing_rate_linear>bins(I))))/(1+0*sqrt(length(firing_rate(find(predicted_firing_rate_linear<=bins(I+1) & predicted_firing_rate_linear>bins(I))))));
end
mean_rate(find(bins<=-40))=0;

p=polyfit(bins(1:end-1),mean_rate,6);
predicted_Dyna_N_Linear=polyval(p,predicted_firing_rate_linear);
pred=polyval(p,bins(1:end-1));
predicted_rate_NO_DNL=polyval(p,transformed_stimulus_linear);
BETA=nlinfit(bins(1:end-1),mean_rate,@sigfun,[mean_rate(1) mean_rate(end) mean(predicted_firing_rate_linear) 1]);
predicted_Dyna_N_Linear=sigfun(BETA,predicted_firing_rate_linear);
predicted_Dyna_N_Linear(find(predicted_Dyna_N_Linear<0))=0;

time1=1;%max(find(time<=55));
time2=length(firing_rate);%max(find(time<=60));

VAF_Dynamicnonlinear=1-var(firing_rate(time1:time2)-predicted_Dyna_N_Linear(time1:time2))/var(firing_rate(time1:time2))
VAF_NO_DNL=1-var(firing_rate(time1:time2)-predicted_rate_NO_DNL(time1:time2))/var(firing_rate(time1:time2))

figure
plot(stimulus(time1:time2)),hold on, plot(firing_rate(time1:time2)), plot(predicted_Dyna_N_Linear(time1:time2))

RR=[];RR=corrcoef(stimulus,predicted_Dyna_N_Linear);
R=[];R=abs(RR(1,2));


