function [ bcn, spike_burst ] = burst_checker( ti, nspikes, dts, dtq )
%burst_checker
%   burst_checker takes inputs
%       ti = vector of length nspikes with all the spike times
%       nspikes = number of spikes input
%       dts = maximum distance in time between successive spikes in burst
%       dtq = minimum quiescent time before first spike of burst
%  returns the following values
%        bcn = number of bursts or burst count
%        spike_burst = vector for every spike which identifies if spike is
%           part of burst and the burst number of each burst



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Parameters of burst_check
%dtq = 100; %delta time for quiescence... in order to be considered burst
%dts = 4;  %time in ms the successive spikes must be less than apart in order to be part of burst

%End parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
%Data Structures

nb = 0;  %total number of bursts
begin_burst  = zeros(nspikes,1);

dt = zeros(nspikes,1);  %time between currents spike and last spike
spike_burst = zeros(nspikes,1);  %0 if not part of burst, j if part of burst j
quiet = zeros(nspikes,1);  %identifies if long enough quiescent period happend
%spike_burst_order = zeros(nspikes,1); %0 if not part of burst, k if part of burst of size k

%End Data Structures
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%calculate dt
dt(1) = dtq + ti(1);  %assuming quiescence from initial
for i = 2:nspikes
    dt(i) = ti(i) - ti(i-1);
end

%calculate if preceding period is quiet

for i = 1:nspikes
    if dt(i) > dtq
        quiet(i) = 1;
    end 
end

%start at each spike that start with quiescent and see if next spike is
%less than dts
bcn = 0;  %burst counter
for i = 1:nspikes - 1
    if quiet(i) == 1
        
        if dt(i + 1) < dts
           bcn = bcn + 1;
           begin_burst(i) = bcn;
           spike_burst(i) = bcn;
           spike_burst(i+1) = bcn;
        end
    end
end

%now go through bursts and see if there are more successive spikes in burst
for i = 1:nspikes-1
    if spike_burst(i) > 0 && dt(i + 1) < dts
        spike_burst(i+1) = spike_burst(i);
    end
end

%
        
        

end %function burst_checker.m

