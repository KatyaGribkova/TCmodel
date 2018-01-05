classdef Neuron_Regular_Spiking < handle
    %Re do to model regular spiking of paper, units will be
    % time is ms
    % potential is in mV
    % conductance is in micro S / cm^2
    % current will be nano amps
    % capacitance will be in nano Farrads
    
    
    properties
        Model = 'Neuron_L4';
%         g_L =  2.89*10^-2;  %loss conductance uS/cm^2
%         E_L = -70;          %loss potential mV
%         g_Na = 1.44*10;     %fast sodium conductance, ms/ mm^2
%         E_Na = 50 ;         %               ms/mm^2
%         g_K = 1.44;       %   potassium rectifier
%         E_K = -90;         %    potassium potential 
%         g_M = 2.0276*10^-2;   %   slow inactivativing potassium
%         tau_p = 1000  ;    %    timescale of slow potassium
        
        g_L =  4.8128;     %loss conductance nS
        E_L = -60.2354;          %loss potential mV
       g_Na = 3000;     %fast sodium conductance, nS
        E_Na = 50 ;         %               ms/mm^2
       g_K = 140;       %   potassium rectifier
        E_K = -90;         %    potassium potential 
        g_M = 1.5;   %   slow inactivativing potassium
        tau_p = 180  ;    %    timescale of slow potassium
        
                     %length of neurong
                 %surface area of neuron


%         cm = .29;        %capacitance of membrane/ surface area
%         Ie_A = -000;     %current injected in nA/mm^2
        
        cm = 109.3865;        %capacitance of pF
        Ie_A = -000;     %current injected in nA/mm^2
        
        NumSyn = 0;  %Number of synapses
        Syn               %Synapse matrix
    

        m = 0.0101;       %activation of fast sodium channel
        h = 0.9659;       %de-inactivation fast sodium channel
        n = 0.1559;       %potassium rectifier channel
        p = 0;            %slow potassium channel
%         
% 
%         htau = 1;
%         mtau = 1;
%         ntau = 1;
%         ptau = 1;
%         jHtau = 1;
%         Vtau = 1;



        V = -60;         %intial potential of membrane (mV)
        %Vx = -10;          % shift of firing potential
        Vt = -47.9;         %shift threshhold of firing


        s = 0.;          %T channel activation
        u = 0.1;            %T channel de-inactivation
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Synapse Input set at 0
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        EsGs = 0;   %Sum of synpase potential and conductance from synpases
        Gs = 0;     %Sum of conductances (g_s)
        
%         mavxe = exp(-(- Vt - 13)/4);
%         mbvxe = exp((- Vt - 40)/5);
% 
%         havxe = exp(-(- Vt - 17)/18);
%         hbvxe = exp(-(- Vt - 40)/5);
% 
%         navxe = exp(-(- Vt - 15)/5);
%         nbvxe = exp(-(- Vt - 10)/40);
% 
%         pivxe = exp(- 35/10);
%         pavxe = exp(35/20);
%         pbvxe = exp(-35/20);


        mavxe = exp(-(47.9 - 13)/4);
        mbvxe = exp((47.9 - 40)/5);

        havxe = exp(-(47.9 - 17)/18);
        hbvxe = exp(-(47.9 - 40)/5);

        navxe = exp(-(47.9 - 15)/5);
        nbvxe = exp(-(47.9 - 10)/40);

        pivxe = exp(- 35/10);
        pavxe = exp(35/20);
        pbvxe = exp(-35/20);
        
        
    end
    
   
   
     
    methods
        function nn = Neuron(input)
            if(nargin > 0)
                nn.Model = input;
            end
        end
        
   
        
        
        function integrate(nn,dt,EsGs,Gs)
            %This will take one time step of dt and update all variables of
            %neuron
            %EsGs is the sum of the conductances and potential
            %Gs is sum of conductances
            %fast sodium channel variable m
            
            ve = exp(nn.V);
            ve5 = ((ve)^(1/5));
            ve20 = ((ve)^(1/20));
            
            

            
            
            alpha = (-0.328*(nn.V - nn.Vt - 13))/(nn.mavxe/((ve)^(1/4)) - 1);
%             alpha = (-0.328*(nn.V - nn.Vt - 13))/...  %openning
%                  (exp(-(nn.V - nn.Vt - 13)/4) - 1);
            beta = (0.28*(nn.V - nn.Vt - 40))/(nn.mbvxe*ve5 - 1);  %clossing
%             beta = (0.28*(nn.V - nn.Vt - 40))/...
%                 (exp((nn.V - nn.Vt - 40)/5) - 1);  %clossing
            mtau = 1./(alpha + beta);
            m_inf = alpha / (alpha + beta);
    
            %newval = timeint(m_inf,nn.m,dt,mtau);
            nn.m = timeint(m_inf,nn.m,dt,mtau);
    
        %   fast sodium channel (variable h)
        

                        
            alpha = 0.128*nn.havxe/((ve)^(1/18));
%             alpha = 0.128*exp(-(nn.V - nn.Vt - 17)/18);
            beta = 4/(1 + nn.hbvxe/ve5);
%             beta = 4/(1 + exp(-(nn.V - nn.Vt - 40)/5));
            htau = 1./(alpha + beta);
            h_inf = alpha / (alpha + beta);
     
            %newval = timeint(h_inf,nn.h,dt,htau);
            nn.h = timeint(h_inf,nn.h,dt,htau);
       
         %potassium rectifier (variable n)
         

            
            
            
            alpha = (-.032*(nn.V - nn.Vt - 15))/(nn.navxe/ve5 - 1);
%             alpha = (-.032*(nn.V - nn.Vt - 15))/...
%                 (exp(-(nn.V - nn.Vt - 15)/5) - 1);\
            beta = 0.5*nn.nbvxe/((ve)^(1/40));
%             beta = 0.5*exp(-(nn.V - nn.Vt - 10)/40);
            ntau = 1./(alpha + beta);
            n_inf = alpha / (alpha + beta);
%     
            %newval = timeint(n_inf,nn.n,dt,ntau);
            nn.n = timeint(n_inf,nn.n,dt,ntau);
     
            %  potassium slow current
            

            
            p_inf = 1/(1 + nn.pivxe/((ve)^(1/10)));
%             p_inf = 1 /(1 + exp(-(nn.V + 35)/10));
            ptau = nn.tau_p/(3.3*nn.pavxe*ve20 + nn.pbvxe/ve20);
%             ptau = nn.tau_p/(3.3 * exp((nn.V + 35)/20) + exp(-(nn.V + 35)/20));
            
            %newval = timeint(p_inf,nn.p, dt, ptau);
            nn.p = timeint(p_inf,nn.p, dt, ptau);
            

         
%     %integrate potential
            gg_Na = nn.g_Na*nn.h*nn.m^3;
            gg_K = nn.g_K*nn.n^4;
            gg_M = nn.g_M*nn.p;
    
          Vinf = (nn.g_L*nn.E_L + gg_Na*nn.E_Na + ...
            gg_K*nn.E_K + gg_M*nn.E_K + nn.Ie_A + EsGs)/...        %Input from synpase
            (nn.g_L + gg_Na + gg_K + gg_M + Gs);
    
          Vtau = nn.cm/...
            (nn.g_L + gg_Na + gg_K + gg_M+ Gs);

            %newval = timeint(Vinf,nn.V,dt,Vtau);
            nn.V = timeint(Vinf,nn.V,dt,Vtau);
%             
            
        end
    end
end
    

