classdef Neuron_FS < handle
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Based of parametersdeleuze2012 paper
    %Re do to model fast spiking of paper, units will be
    % time is ms
    % potential is in mV
    % conductance is in nS
    % capacitance will be in pFarrads
    % current = nA
   
    
    properties
        Model = 'Fast Spiking';
         g_L =  3.9;  %loss conductance nS
         E_L = -58.0;          %loss potential mV
        g_Na = 1500;     %fast sodium conductance, ms
        E_Na = 50 ;         %               ms/mm^2
         g_K = 440;      %   potassium rectifier
         E_K = -100;          %    potassium potential 
         g_H =1.3445;  %matching h-current to leak
         E_H = -33;              %potential of h-current
      
                 
        cm = 70.;        %capacitance of membrane
        Ie_A = -000;     %current injected in nA/mm^2
        
       
        m = 0.0126;       %activation of fast sodium channel
        h = 0.9964;       %de-inactivation fast sodium channel
        n = 0.0309;       %potassium rectifier channel
%         

        V = -58.66;         %intial potential of membrane (mV)
        Vx = -32.;          % shift of firing potential

        jH = 0.2950;
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Synapse Input set at 0
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        EsGs = 0;   %Sum of synpase potential and conductance from synpases
        Gs = 0;     %Sum of conductances (g_s)
        
        
        
    end
    
   
    methods
        function nn = Neuron(input)
            if(nargin > 0)
                nn.Model = input;
            end
        end
        

        function integrate(nn,dt,EsGs, Gs)    
             alpha = 1.3*0.32*(13 - nn.V + nn.Vx)/...
                 (exp((13- nn.V + nn.Vx)/4)-1);%openning
           
            beta = 1.4*0.28*(nn.V - nn.Vx -40)/...
                (exp((nn.V - nn.Vx -40)/5)-1);  %clossing
            mtau = 1./(alpha + beta);
            m_inf = alpha / (alpha + beta);
%     
            newval = timeint(m_inf,nn.m,dt,mtau);
             nn.m = newval;
%     
        %   fast sodium channel (variable h)
            alpha = 1.3*0.128*exp((17 -nn.V + nn.Vx)/18);
            beta = 1.3*4/(1 + exp((40 - nn.V + nn.Vx)/5));
            htau = 1./(alpha + beta);
            h_inf = alpha / (alpha + beta);
     
            newval = timeint(h_inf,nn.h,dt,htau);
            nn.h = newval;
%        
%          %potassium rectifier (variable n)
            alpha = 1.4*0.032*(15 - nn.V+ nn.Vx)/(exp((15 - nn.V+ nn.Vx)/5)-1);
            beta = 1.6*0.5*exp((10 - nn.V + nn.Vx)/40);
            ntau = 1./(alpha + beta);
            n_inf = alpha / (alpha + beta);
%     
            newval = timeint(n_inf,nn.n,dt,ntau);
            nn.n = newval;
     

             %H-current model (originally from McCormick and Hugenard)
             jH_inf = 1./(1 + exp((nn.V + 75)/(5.5)));
             jHtau = 1/(exp(-14.59 - 0.086*nn.V) + ...
                 exp(-1.87 + 0.0701*nn.V));
%             
             newval = timeint(jH_inf, nn.jH, dt, jHtau);
             nn.jH = newval;
%            

         
%     %integrate potential
             gg_Na = nn.g_Na*nn.h*nn.m^3;
             gg_K = nn.g_K*nn.n^4;
%              gg_T = nn.g_T*nn.H*nn.M^2;
             gg_H = nn.g_H*nn.jH;
    
          Vinf = (nn.g_L*nn.E_L + gg_Na*nn.E_Na + ...
            gg_K*nn.E_K  + gg_H*nn.E_H + ...
            nn.Ie_A + EsGs)/...        %Input from synpase
            (nn.g_L + gg_Na + gg_K + gg_H + Gs);
    
          Vtau = nn.cm/...
            (nn.g_L + gg_Na + gg_K +gg_H + Gs);

            newval = timeint(Vinf,nn.V,dt,Vtau);
            nn.V = newval;
%             
            
        end
    end
end
    

