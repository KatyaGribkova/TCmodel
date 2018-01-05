%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Include recipricol loop in TRN simulation:  see notes from Feb 3
%Uses Class Neuron_Thalamic_Relay
%           Neuron_Fast_Spike
%           Synapse
%           Synapse_Depression
% and uses functions poisson.m
%                    spike_check.m
%                    timeint.m
%                    burst_checker.m
%
%2 neuron case:  
%                
%             TC:  Neuron1 thalamocortical cell linked to Cortex and TRN,
%                  also recieves inputs from Input1 (excite) and Input2
%                  (inhib)
%             TRN: Neuron 2 reticular neuron with excitatory input from TC and sends inhibiator input to
%                   TC 
%                   and inhibitory input from TRN (input 1), using saturating synapse.
%             Cortex: Neuron 3 Corticol (Fast spiking) neuron with depressing
%               (excitatory)  synapse inputs from Neuron3.
%                
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
%Two neuron case, with A exciting B and B inhibiting C.  Both A and C will
%also recieve outside input.


NN = 3;                     % Number of neurons

SimTime =100000;              %Simulation time in ms 
dt = .1 ;                   %time step in ms
Nt = SimTime/dt;        %Number of steps in simulation.

SpkTh = 00;              %Spike Threshold

%SimNum = 1;              %Number of simulations per input frequency

maxrateI = .1;               %max rate of input for inhibit in ms
maxrateE = .1;               %max rate of input for excite in ms

buffer = 2;             %ratio of size of buffer to number of spikes

SimNum = 1;

Vinit = -66;                %setting initial value of each neuron

EqlbTime = 200;                  %time in ms to equalibrate and not gather statistics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Define synapses in matrix Syn.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Excitatory external synapse
SynInE = Synapse_Depression;   %  external excitation synapase
SynInE.E_s = 0;
SynInE.gs = 40;
SynInE.gs = 32;
SynInE.U_se = 0.76;
SynInE.tau_inact = 2.64;
SynInE.tau_rec = 250.0;
SynInE.tau_rec = 125.0;
%SynInE.tau_rec = 20.0;

%Excitatory external synapse for TRN
SynInETRN = Synapse_Depression;
SynInETRN.E_s = 0;
SynInETRN.gs = 16;
SynInETRN.gs = 32;
SynInETRN.U_se = 0.3;
SynInETRN.tau_inact = 10.58;
SynInETRN.tau_rec = 40;  %fixed value for excitation


%Inhibatory synapse from TRN to TC 
SynTRN_TC = Synapse_Depression;
SynTRN_TC.E_s = -100;
SynTRN_TC.E_s = -80;
SynTRN_TC.gs =  20; 
SynTRN_TC.gs =  80;%%%%%%%%%%%%%adjustable parameter
SynTRN_TC.U_se = 0.62;
SynTRN_TC.tau_inact = 16.62;
%SynTRN_TC.tau_inact = 24.62;
SynTRN_TC.tau_rec = 167.29;


% %Excitatory Synapes from TC to TRN
% SynTC_TRN = Synapse_Depression;
% SynTC_TRN.E_s = 0;
% SynTC_TRN.gs =  2*40;
% SynTC_TRN.U_se = 0.76;
% SynTC_TRN.tau_inact = 2.64;
% SynTC_TRN.tau_rec = 562.13;

%Excitatory Synapse (Depressing) TC to C
SynTC_C = Synapse_Depression;
SynTC_C.E_s = 0;
SynTC_C.gs = 20;
SynTC_C.gs = 6675;
SynTC_C.gs = 50;
SynTC_C.U_se = 0.8113;
SynTC_C.U_se = 0.8113;
SynTC_C.tau_inact = 11.52;
SynTC_C.tau_inact = 3;
SynTC_C.tau_rec = 160 ;%1000
SynTC_C.tau_rec = 165;


%TC to FS excitatory synapse (currently copyin SynTC_C)
SynTC_FS = Synapse_Depression;
SynTC_FS.E_s = 0;
SynTC_FS.gs = 70;
SynTC_FS.gs = 0;
SynTC_FS.U_se = 0.28;
SynTC_FS.tau_inact = 5.652;
SynTC_FS.tau_rec = 984 ;%1000

%FS to C inhibitory synapse  (currently copying SynTRN_TC)
SynFS_C = Synapse_Depression;
SynFS_C.E_s = -100;
SynFS_C.gs =  30;             %%%%%%%%%%%%%adjustable parameter
SynFS_C.gs =  0;
SynFS_C.U_se = 0.2;
SynFS_C.tau_inact = 7.162;
SynFS_C.tau_rec = 511.41;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Create Neurons
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

TC = Neuron_Thalamic_Relay_deleuze2012; %Thalamocortical cell
TC.g_T=40;

TRN = Neuron_TRN;  %Thalamic Reticular Nucleaus

FS = Neuron_FS;   %Fast spiking interneuron


C = Neuron_Regular_Spiking;  %Cortex
%C.g_M=1.5;
%C.tau_p=9000;
%C.E_L= -80

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %Allocate space for keeping track of percent transmission
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Spike_E_his = zeros(SimNum,SimNum);
Spike_I_his = zeros(SimNum,SimNum);
Spike_TC_his = zeros(SimNum,SimNum);
Spike_C_his = zeros(SimNum,SimNum);
BurstTC_his = zeros(SimNum,SimNum);

Spike_E_Times = zeros(SimNum,SimNum,round(maxrateE*SimTime*buffer));
Spike_I_Times = zeros(SimNum,SimNum,round(maxrateE*SimTime*buffer));

Spike_TC_Times = zeros(SimNum,SimNum,round(maxrateE*SimTime*buffer));
Spike_C_Times = zeros(SimNum,SimNum,round(maxrateE*SimTime*buffer));
Spike_TRN_Times = zeros(SimNum,SimNum,round(maxrateE*SimTime*buffer));

% V_E_Traces = zeros(SimNum,SimNum,Nt);
% V_I_Traces = zeros(SimNum,SimNum,Nt);
% V_TC_Traces = zeros(SimNum,SimNum,Nt);
% V_TRN_Traces = zeros(SimNum,SimNum,Nt);
% V_C_Traces = zeros(SimNum,SimNum,Nt);

MI_E_TC = zeros(SimNum,SimNum);
MI_TC_C = zeros(SimNum,SimNum);
MI_E_C = zeros(SimNum,SimNum);
%MI_I_TRN = zeros(SimNum,SimNum);
%MI_TRN_TC = zeros(SimNum,SimNum);
MI_TRN_C = zeros(SimNum,SimNum);



% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %Allocate Quantities for optimization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
V_C = zeros(Nt,1);     %Voltage of C, TC, TRN for debugging
V_C(1) = Vinit;
V_TC = zeros(Nt,1);
V_TC(1) = Vinit;
V_TRN = zeros(Nt,1);
V_TRN(1) = Vinit;
V_FS = zeros(Nt,1);
V_FS(1) = Vinit;
MIcount=0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Debugging storage
%     ESynInI = zeros(Nt,1);
%     ESynInE = zeros(Nt,1);
%     ESynTC_TRN = zeros(Nt,1);
%     ESynTRN_TC = zeros(Nt,1);
%     ESynTC_C = zeros(Nt,1);
% 
%    
%END Debugging
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
inputE = zeros(Nt,1);%Modeled spike excitatory
inputI = zeros(Nt,1);%Modeled spike INhibitory
time = zeros(Nt,1);   %time vector for plotting
% 
% rate1His = zeros(SimNum,1);    %history of Isi2 and
% rate3His = zeros(SimNum,1);
% 
% TransR(:,:) = 0;
% TransN1N3(:,:) = 0;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %Generate Input
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 


run=0;
%TRNrates=[0 0.00063 0.00078 0.00098 0.00122 0.00153 0.00191 0.00238 0.00298 0.00373 0.00466 0.00582 0.00728 0.00909 0.01137 0.01421 0.01776 0.02220 0.02776 0.03469 0.04337 0.05421 0.06776 0.08470 0.10588 0.13235 0.16544 0.20680];
%TCrates=[0.00050 0.00063 0.00078 0.00098 0.00122 0.00153 0.00191 0.00238 0.00298 0.00373 0.00466 0.00582 0.00728 0.00909 0.01137 0.01421 0.01776 0.02220 0.02776 0.03469 0.04337 0.05421 0.06776 0.08470 0.10588 0.13235 0.16544 0.20680];
TRNrates=[0.00063];
TCrates=[0.00122];
 for l = 1:SimNum     %l varies rate3
    for m = 1:SimNum    %m vaires rate1
        %l*m
        %run=run+1
%          rate1 = (m-1)*(maxrateI/(SimNum-1))
         %rate1 = (m-1)*(maxrateI/(SimNum-1));
         %rate3 = (l)*(maxrateE/(SimNum));
         rate1 = TRNrates(m);
         rate3 = TCrates(l);
         %rate3His(l) = rate3;
         %rate1His(m) = rate1;% %         Isi1 = round(minIsi1 + (maxIsi1/SimNum)*(m-1))
% %         Isi2 = round(minIsi2 + (maxIsi2/SimNum)*(l-1))
% %         Isi2His(l) = Isi2;
% %         Isi1His(m) = Isi1;
        %spikein =0;             %number of spikes received to  neuron 3
        %spikeTC = 0;           %number of spikes of neuron 3
        n1spikeout = 0;         %number of spikes output of neuron1
%     
     inputE(:) = -60;
     inputI(:) = -60;

     %reset potential
   
    TC.V = Vinit;
    C.V = Vinit;
    TRN.V = Vinit;
    FS.V = Vinit;
%    
%    % reset synpases to closed
   
    SynInETRN.E = 0.0;
    SynInETRN.R = 1.0;
    SynInETRN.I = 0;
   
    SynInE.E = 0.0;
    SynInE.R = 1.0;
    SynInE.I = 0;
    
    SynTC_C.E = 0.0;
    SynTC_C.R = 1.0;
    SynTC_C.I = 0;
    
    SynTRN_TC.E = 0.0;
    SynTRN_TC.R = 1.0;
    SynTRN_TC.I = 0;
    
    SynTC_FS.E = 0.0;
    SynTC_FS.R = 1.0;
    SynTC_FS.I = 0.;
    
    SynFS_C.E = 0.0;
    SynFS_C.R = 1.;
    SynFS_C.I = 0;
    
    
%     SynTC_TRN.E = 0.0;
%     SynTC_TRN.R = 1.0;
%     SynTC_TRN.I = 0;
    

   
    %Input spike generation

%     
% %     %generate spike train of rate1 with function poisson.m
     if rate1 > 0
         [nspikeI,tempspikeI] = poisson(rate1,EqlbTime,SimTime);  
       %  interval=100;tempspikeI=500:interval:2500-interval;nspikeI=length(tempspikeI);
     else    nspikeI = 0;
     end
    if rate3 > 0
         [nspikeE,tempspikeE] = poisson(rate3,EqlbTime,SimTime);
        % interval=50;tempspikeE=500:interval:2500-interval;nspikeE=length(tempspikeE);
    end
     spkIcnt = 1;   %spike counter I(looking for first spike)
     spkEcnt = 1;   %spike counter E
     spikein = 0;   %input spike number
     
     spikeC = 0;  %cortex spike counter
     spikeTC = 0;
     spikeTRN = 0;
     spikeE = 0;
     spikeI = 0;
%     
%     %reset inputs to -60 mV
     inputE(:) = -60;
     inputI(:) = -60;
%     
     %go through spike list and place a spike at t(i+1) if spike is found
     %between times t(i) and t(i+1)
     for i = 1:Nt 
         time(i) = (i-1)*dt;
         if spkIcnt <= nspikeI   %check bounds of tempspike1
             if tempspikeI(spkIcnt) < (i+1)*dt
                 inputI(i+1) = 30;   
                 spkIcnt = spkIcnt + 1;
             end
         end
%         
         if spkEcnt <= nspikeE   %check bounds of tempspike3
             if tempspikeE(spkEcnt) < (i+1)*dt
                 inputE(i+1) = 30;
                 spkEcnt = spkEcnt + 1;
                 spikein = spikein + 1;
             end
         end
     end
% %

%inputE(:)=-60; %used to zero out input
     clear tempspike1; clear tempspike2  %clear memory of tempspikes    
% % 
% % 
% %  
% % 

fprintf('Integrating... \n')

for i = 2:Nt
    
     %Integrate
    
     EsGsTC = ...
         SynInE.Gs*SynInE.E_s + ...   %excit input to TC
         SynTRN_TC.Gs*SynTRN_TC.E_s;  %feedback TRN input

     EsGsC = SynTC_C.Gs*SynTC_C.E_s + ...      %TC input to C
            SynFS_C.Gs*SynFS_C.E_s;
      
        EsGsTRN = SynInETRN.Gs*SynInETRN.E_s ;   %input to TRN
        
        EsGsFS = SynTC_FS.Gs*SynTC_FS.E_s;      %input to FS

%     
    
     GsTC = SynInE.Gs  + SynTRN_TC.Gs ;
     GsC = SynTC_C.Gs + SynFS_C.Gs;
     GsTRN = SynInETRN.Gs;
     GsFS = SynTC_FS.Gs;
%     
     
     TC.integrate(dt,EsGsTC,GsTC);
     TRN.integrate(dt,EsGsTRN,GsTRN);
     C.integrate(dt,EsGsC,GsC);
     FS.integrate(dt,EsGsFS,GsFS);
    
     %integrate synapses
    SynInE.integrate(dt);
    SynTRN_TC.integrate(dt);
    SynInETRN.integrate(dt);
    SynTC_C.integrate(dt);
    SynFS_C.integrate(dt);
    SynTC_FS.integrate(dt);
     
    
%     
     V_TC(i) = TC.V;
     V_TRN(i) = TRN.V;
     V_C(i) = C.V;
     V_FS(i) = FS.V;
     
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%      %Debugging Synapses
%    
%     ESynInI = SynInI.E;
%     ESynInE = SynInE.E;
%     ESynTC_TRN = SynTC_TRN.E;
%     ESynTRN_TC = SynTRN_TC.E;
%     ESynTC_C = SynTC_C.E;
    %END Debugging Synpases
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
     %process spiking (delayed one dt in time)
     %check only after time > 3 
     delay=10;
      if i > 10 %was 3  Needs to be greater than delay or else it will look for negative time
     if spike_check(inputI(i-delay:i),SpkTh)  %External excitation of TRN
         SynInETRN.SpikeReceived;
         spikeI = spikeI + 1;
         Spike_I_Times(l,m,spikeI) = i*dt;
     end
     
     if spike_check(inputE(i-delay:i),SpkTh)  %excite input to TC
         SynInE.SpikeReceived;
         spikeE = spikeE + 1;
         Spike_E_Times(l,m,spikeE) = i*dt;
         
     end
%      
     if spike_check(V_TC(i-delay:i),SpkTh)    %keeps track of spikes of TC to C
         SynTC_C.SpikeReceived;
         SynTC_FS.SpikeReceived;
         spikeTC = spikeTC + 1;
         Spike_TC_Times(l,m,spikeTC) = i*dt;
     end
% 
    if spike_check(V_TRN(i-delay:i),SpkTh)   %spike from TRN 
        SynTRN_TC.SpikeReceived;
        spikeTRN = spikeTRN + 1;
        Spike_TRN_Times(l,m,spikeTRN) = i*dt;
    end
%     
    if spike_check(V_C(i-delay:i),SpkTh)   %spike from C
        spikeC = spikeC +1;
        Spike_C_Times(l,m,spikeC) = i*dt;
    end
    
    if spike_check(V_FS(i-delay:i),SpkTh)  %spike from FS
        SynFS_C.SpikeReceived;
    end
%      
      end
    
%%%%%%%%%%%%%%Do Burst Checking%%%%%%%%%%%%%%%%%%%%%



 end %time stepping
 
%  V_E_Traces(l,m,:) = InputE(:);
%  V_I_Traces(l,m,:) = InputI(:);
%  V_TC_Traces(l,m,:) = V_TC(:);
%  V_TRN_Traces(l,m,:) = V_TRN(:);
%  V_C_Traces(l,m,:) = V_C(:);
 

Spike_E_his(l,m) = spikein;
Spike_TC_his(l,m) = spikeTC;
Spike_C_his(l,m) = spikeC;
Spike_I_his(l,m) = nspikeI;


% MI_E_TC(l,m) = AIMIE(Spike_E_Times(l,m,:),Spike_TC_Times(l,m,:));
% MI_TC_C(l,m) = AIMIE(Spike_TC_Times(l,m,:),Spike_C_Times(l,m,:));
% MI_E_C(l,m) = AIMIE(Spike_E_Times(l,m,:),Spike_C_Times(l,m,:));
% MI_TRN_C(l,m) = AIMIE(Spike_TRN_Times(l,m,:),Spike_C_Times(l,m,:));

%save(mfilename)

%MIcount=MIcount+1

    end
 end
% % 
% % TransR(l,m) = spikeTC / spikein;    %keeping track of transmission rate
% % TransN1N3(l,m) = spikeTC/n1spikeout; %keeps track of spike ratio
% % 
% % end   %varying Isi1
% % end   % varying Isi2
% % 
% % % plot(TransR)
% % % title('Transmission Rate %')

% 

 timecourse=((dt:dt:SimTime)/1000);
 subplot(5,1,1);plot(timecourse,V_C)
 title('cortex')
 subplot(5,1,2);plot(timecourse,V_TC)
 title('TC')
 subplot(5,1,3);plot(timecourse, inputE)
 title('input to TC')
 subplot(5,1,4);plot(timecourse,V_TRN)
 title('TRN')
 subplot(5,1,5);plot(timecourse,inputI)
 title('Input to TRN')

 
% 

