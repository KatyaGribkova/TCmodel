classdef Neuron_TRN < handle
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Based of parametersdeleuze2012 paper
    %Re do to model fast spiking of paper, units will be
    % time is ms
    % potential is in mV
    % conductance is in nS
    % capacitance will be in pFarrads
    % current = nA
   
    
    properties
        Model = 'TRN';
         g_L =  3.7928;  %loss conductance nS
         E_L = -57.0;          %loss potential mV
        g_Na = 3000;     %fast sodium conductance, ms
        E_Na = 50 ;         %               ms/mm^2
         g_K = 400;      %   potassium rectifier
         E_K = -100;          %    potassium potential 
        g_T = 21.  ;  % t-curret calcium
        E_T = 106.7 ;           %potential of t current
         g_H =.0192;  %matching h-current to leak
         E_H = -33;              %potential of h-current
         
        
         g_KS = 3.5;      %conductance of non-inactivating current slow potassium
         tau_KS = 200;  %time constant
      
                     %length of neurong
                 %surface area of neuron

        Q = 2;    %variable for temperature and Ca++ concentration (Coulter 1989)

        cm = 75.0;        %capacitance of membrane
        Ie_A = -000;     %current injected in nA/mm^2
        
       
        m = 0.0126;       %activation of fast sodium channel
        h = 0.9964;       %de-inactivation fast sodium channel
        n = 0.0309;       %potassium rectifier channel
%         


%         htau = 1;
%         mtau = 1;
%         ntau = 1;
%         ptau = 1;
%         Mtau = 1;
%         Htau = 1;
%         jHtau = 1;
%         Vtau = 1;

        V = -65.4761;         %intial potential of membrane (mV)
        Vx = -48.;          % shift of firing potential
        Vt =8;         %shift threshhold for t-current

nu
%         s = 0.;          %T channel activation
%         u = 0.1;            %T channel de-inactivation

        M = .2044;             % activation of T-current
        H = .0198 ;            % inactivation
        
        % non-inactivating slow potassium
        KS = .01;              %slow potassium activation
        %H current addition
        jH = 0.2950;
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Synapse Input set at 0
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        EsGs = 0;   %Sum of synpase potential and conductance from synpases
        Gs = 0;     %Sum of conductances (g_s)
        
        
%         mavxe = exp((13 + Vx)/4);
%         mbvxe = exp((Vx + 40)/5);
%         
%         havxe = exp((17 + Vx)/18);
%         hbvxe = exp((40 + Vx)/5);
%         
%         navxe = exp((15 + Vx)/5);
%         nbvxe = exp((10 + Vx)/40);        
%         
%         Mivxe = exp(-(- Vt + 57)/6.2);
%         Mavxe = exp(-(-Vt +132)/16.7);
%         Mbvxe = exp((- Vt + 16.8)/18.2);
%         
%         Hivxe = exp((- Vt + 81)/4);
%         Hlessvxe = exp(( - Vt + 467)/66.6);
%         Hmorevxe = exp(-(- Vt + 22)/10.5);
% 
%         jHivxe = exp((75)/(5.5));
%         jHavxe = exp(-14.59);
%         jHbvxe = exp(-1.87);
%         
%         pivxe = exp(- 35/10);
%         pavxe = exp(35/20);
%         pbvxe = exp(-35/20);

        mavxe = exp((13 - 48.)/4);
        mbvxe = exp((-48. + 40)/5);
        
        havxe = exp((17 - 48.)/18);
        hbvxe = exp((40 - 48.)/5);
        
        navxe = exp((15 - 48.)/5);
        nbvxe = exp((10 - 48.)/40);        
        
        Mivxe = exp(-(- 8 + 57)/6.2);
        Mavxe = exp(-(-8 +132)/16.7);
        Mbvxe = exp((- 8 + 16.8)/18.2);
        
        Hivxe = exp((- 8 + 81)/4);
        Hlessvxe = exp(( - 8 + 467)/66.6);
        Hmorevxe = exp(-(- 8 + 22)/10.5);

        jHivxe = exp((75)/(5.5));
        jHavxe = exp(-14.59);
        jHbvxe = exp(-1.87);
        
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
        

        function integrate(nn,dt,EsGs, Gs)   
            
            
             ve = exp(nn.V);
             ve5 = ((ve)^(1/5));
             ve4 = ((ve)^(1/4));
             ve20 = ((ve)^(1/20));
 
             
            alpha = 1.3*0.32*(13 - nn.V + nn.Vx)/((nn.mavxe/ve4) - 1);            
%             alpha = 1.3*0.32*(13 - nn.V + nn.Vx)/...
%                  (exp((13- nn.V + nn.Vx)/4)-1);%openning
           
            beta = 1.4*0.28*(nn.V - nn.Vx -40)/((ve5 / nn.mbvxe)-1);
%             beta = 1.4*0.28*(nn.V - nn.Vx -40)/...
%                 (exp((nn.V - nn.Vx -40)/5)-1);  %clossing
            mtau = 1./(alpha + beta);
            m_inf = alpha / (alpha + beta);
%     
            %newval = timeint(m_inf,nn.m,dt,mtau);
             nn.m = timeint(m_inf,nn.m,dt,mtau);
%     
        %   fast sodium channel (variable h)
        
        
        
            alpha = 1.3*0.128*(nn.havxe/((ve)^(1/18)));        
%             alpha = 1.3*0.128*exp((17 -nn.V + nn.Vx)/18);
            beta = 1.3*4/(1 + (nn.hbvxe/ve5));
%             beta = 1.3*4/(1 + exp((40 - nn.V + nn.Vx)/5));
            htau = 1./(alpha + beta);
            h_inf = alpha / (alpha + beta);
     
            %newval = timeint(h_inf,nn.h,dt,htau);
            nn.h = timeint(h_inf,nn.h,dt,htau);
%        
%          %potassium rectifier (variable n)



            alpha = 1.4*0.032*(15 - nn.V+ nn.Vx)/((nn.navxe/ve5)-1);
%             alpha = 1.4*0.032*(15 - nn.V+ nn.Vx)/(exp((15 - nn.V+ nn.Vx)/5)-1);
            beta = 1.6*0.5*(nn.nbvxe/((ve)^(1/40)));
%             beta = 1.6*0.5*exp((10 - nn.V + nn.Vx)/40);
            ntau = 1./(alpha + beta);
            n_inf = alpha / (alpha + beta);
%     
            %newval = timeint(n_inf,nn.n,dt,ntau);
            nn.n = timeint(n_inf,nn.n,dt,ntau);
     
            
            %Calcium T-current activation
            
            
            
            
            M_inf = 1 / (1 + nn.Mivxe/((ve)^(1/6.2)));            
%             M_inf = 1 / (1 + exp(-(nn.V - nn.Vt + 57)/6.2));
            Mtau = (0.612 + 1/(nn.Mavxe/((ve)^(1/16.7)) + nn.Mbvxe*((ve)^(1/18.2))))/nn.Q;
%             Mtau = (0.612 + 1/(exp(-(nn.V-nn.Vt +132)/16.7) ...
%                 + exp((nn.V - nn.Vt + 16.8)/18.2)))/nn.Q;
            
            %newval = timeint(M_inf,nn.M,dt,Mtau);
            nn.M = timeint(M_inf,nn.M,dt,Mtau);
            
            %Calcaium T-current inactivation
            
            
            
            H_inf = 1./(1 + nn.Hivxe*ve4);           
%             H_inf = 1./(1 + exp((nn.V - nn.Vt + 81)/4));

            if nn.V < -81
                
               Htau = (nn.Hlessvxe*((ve)^(1/66.6)))/nn.Q;               
%                Htau = (exp((nn.V - nn.Vt + 467)/66.6))/nn.Q;
              
            else
                Htau = (28 + nn.Hmorevxe/((ve)^(1/10.5)))/nn.Q;                
%                 Htau = (28 + exp(-(nn.V - nn.Vt + 22)/10.5))/nn.Q;
          
 
            end
            
            %newval = timeint(H_inf,nn.H,dt,Htau);
            nn.H = timeint(H_inf,nn.H,dt,Htau);

             %H-current model (originally from McCormick and Hugenard)
             
             
             jH_inf = 1./(1 + nn.jHivxe*((ve)^(1/5.5)));             
%              jH_inf = 1./(1 + exp((nn.V + 75)/(5.5)));
             jHtau = 1/(nn.jHavxe/((ve)^(0.086)) + nn.jHbvxe*((ve)^(0.0701)));
%              jHtau = 1/(exp(-14.59 - 0.086*nn.V) + ...
%                  exp(-1.87 + 0.0701*nn.V));
%             
             %newval = timeint(jH_inf, nn.jH, dt, jHtau);
             nn.jH = timeint(jH_inf, nn.jH, dt, jHtau);
             
             %slow non-inactivating potassium current

             
             
             p_inf = 1/(1 + nn.pivxe/((ve)^(1/10)));
%              p_inf = 1/(1 + exp(-(nn.V + 35)/10));
             ptau = nn.tau_KS/(3.3*nn.pavxe*ve20 + nn.pbvxe/ve20);
%              ptau = nn.tau_KS/(3.3*exp((nn.V + 35)/20) + exp(-(nn.V +35)/20));
%            
            %newval = timeint(p_inf,nn.KS,dt,ptau);
            nn.KS = timeint(p_inf,nn.KS,dt,ptau);
         
%     %integrate potential
             gg_Na = nn.g_Na*nn.h*nn.m^3;
             gg_K = nn.g_K*nn.n^4;
             gg_T = nn.g_T*nn.H*nn.M^2;
             gg_H = nn.g_H*nn.jH;
             gg_KS = nn.g_KS*nn.KS;
    
          Vinf = (nn.g_L*nn.E_L + gg_Na*nn.E_Na + ...
            gg_K*nn.E_K + gg_T*nn.E_T + gg_H*nn.E_H + ...
            gg_KS*nn.E_K+...
            nn.Ie_A + EsGs)/...        %Input from synpase
            (nn.g_L + gg_Na + gg_K + gg_T + gg_H + gg_KS+ Gs);
    
          Vtau = nn.cm/...
            (nn.g_L + gg_Na + gg_K + gg_T +gg_H + gg_KS+ Gs);

            %newval = timeint(Vinf,nn.V,dt,Vtau);
            nn.V = timeint(Vinf,nn.V,dt,Vtau);
%             
            
        end
    end
end
    

