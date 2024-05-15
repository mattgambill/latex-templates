function [Range, Velocity, Azimuth] = radarSimulator(R0,theta,v,Np,fp,CPI)
% [matcharrayTx1, matcharrayTx2, rangedopplerTx1, rangedopplerTx2] = radarSimulator(R0,theta,v,Np)
% 
% Simulates Radar for a moving target at velocity v, inital range R0, and initial azimuth theta.
% Uses Np pulses in Coherent Processing Interval CPI.
% Contributors: Matt Gambill, Dr. Ertin

    % Physical Constants
    c = physconst('LightSpeed');            % speed of light

    % Parameters for MISO Radar System
    fc = 10.5e9;        % center frequency = 10.5 GHz
    lambda = c/fc;      % wavelength of radar system
    BW = 125e6;         % total system bandwidth = 125 MHz
    %fp = 1e3;           % pulse repetition frequency = 1 kHz
    Tp = 1/fp;          % pulse repetition interval = 1 ms
    %Np = 64;            % number of pulses in CPI
    Pt = 0.1;            % Watts Transmitted Power
    F = 8;              % db Noise Figure
    STEMP = systemp(F); % System temperature
    beta = 120e6;          % sweep bandwidth = 120 MHz
    tau = 80e-6;        % pulse width = 80 usec
    kb = 1.38064852e-23; % boltzman constant
    G = db2pow(3);      % Antenna Gain
    %BeamWidth = 60*pi/180; % radians
    %M= 64;
    %Tc = Tp*M;             % CPI


    % Parameters for Sampled System
    fs = 250e6;         % sample rate = 250 Msamples/second
    Ts = 1/fs;          % sample period


    % Time Vector
    t = 0:Ts:Np*Tp-Ts;     % time vector (Np+1) Tp long

    % Received signals
    za1 = zeros(size(t));   % received signal from transmitter 1
    za2 = zeros(size(t));   % received signal from transmitter 2


    % Simulate noise
    sigma_n = sqrt(kb*STEMP*BW); % You should calculate using k T0 BW F
    noise = sigma_n*(randn(1,length(t))) + 1i*(sigma_n*(randn(1,length(t))));

    % Radar Postitions
    tx1=[ 0, lambda/4]; tx2=[0, -lambda/4]; rx=[0,0];
    tgt = zeros(Np,2);
    % Calculate azimuth angle, range, and received signals at each location
    for k=0:Np-1
        target=[R0*cos(theta)+k*Tp*v(1) R0*sin(theta)+k*Tp*v(2)];
        tgt(k+1,:) = target;
        Rup1 = norm(tx1-target);
        Rup2= norm(tx2-target);
        Rdown = norm(rx-target);
        % Rup1 = distance between TX1 antenna and target
        % Rup2 = distance between TX2 antenna and target
        % Rdown = distance between RX antenna and target

       Ac1= sqrt(Pt*G^2*lambda^2/((4*pi)^3 * Rup1^2 * Rdown^2));


       Ac2= sqrt(Pt*G^2*lambda^2/((4*pi)^3 * Rup2^2 * Rdown^2));

       za1 = za1 + Ac1*(rpulse(t-k*Tp-(Rup1+Rdown)/c,tau)) .* ...
            (exp(-1i*(2*pi/lambda)*(Rup1+Rdown))) .* ...
            (exp(1i*pi*(beta/tau).*(t-(tau/2)-(k*Tp+(Rup1+Rdown)/c )) .^2));

        za2 = za2 + Ac2*(rpulse(t-k*Tp-(Rup2+Rdown)/c,tau)) .* ...
            (exp(-1i*(2*pi/lambda)*(Rup2+Rdown))) .* ...
            (exp(1i*pi*(beta/tau).*(t-(tau/2)-(k*Tp+(Rup2+Rdown)/c )) .^2));
        % za1 = received signal from TX1
        % za2 = received signal from TX2
        % za = total received signal
    end
    za1n = za1+ noise/2;
    za2n = za2+ noise/2;
    %za = za1n+za2n;
    clear za1 za2 noise
    %Prepare match filter for Tx1
    tm=0:Ts:tau-Ts;
    h1=exp(1i*pi*(beta/tau).*(tm-(tau/2)) .^2); 
    h1=conj(fliplr(h1)); %conjugate and flip the time.
    zeroPad = 1000;
    %
    receivearrayTx1 = zeros(Np,length(tm)+zeroPad);
    receivearrayTx2 = zeros(Np,length(tm)+zeroPad);
    for k=0:Np-1
        receivearrayTx1(k+1,:) = za1n(k*uint64(Tp*fs)+1:round(k*uint64(Tp*fs)+tau*fs+zeroPad));
        receivearrayTx2(k+1,:) = za2n(k*uint64(Tp*fs)+1:round(k*uint64(Tp*fs)+tau*fs+zeroPad));
    end
    %
    clear za1n za2n
    matcharrayTx1=zeros(Np,zeroPad+1);
    matcharrayTx2=zeros(Np,zeroPad+1);
    for k=0:Np-1
        matcharrayTx1(k+1,:) = conv(receivearrayTx1(k+1,:),h1,'valid');
        matcharrayTx2(k+1,:) = conv(receivearrayTx2(k+1,:),h1,'valid');
    end
    %
    clear receivearrayTx1 receivearrayTx2
   

    rangedopplerTx1 = fftshift( fft(matcharrayTx1,[],1),1);
    rangedopplerTx2 = fftshift( fft(matcharrayTx2,[],1),1);
    
    [~, index] = max(abs(rangedopplerTx1), [], 'all','linear');
    [rowtx1, coltx1] = ind2sub(size(rangedopplerTx1),index);
    
    nugrid_norm = 1/Np*(-Np/2:1:(Np/2)-1);
    
    
    tdtx1 = Ts * coltx1;
    % fd_norm = fd/fp
    fdtx1 = nugrid_norm(rowtx1)/fp;
    
    [~, index] = max(abs(rangedopplerTx2), [], 'all','linear');
    [rowtx2, coltx2] = ind2sub(size(rangedopplerTx2),index);
    tdtx2 = Ts * coltx2;
    fdtx2 = nugrid_norm(rowtx2)/fp;

    [Range, Azimuth] = getPos(tdtx1, tdtx2,lambda/2);
    
    Velocity = mean([fdtx1 fdtx2])*lambda/2*fp;
    
    
   %% Plots
   nugrid = 1/Np*(-Np/2:1:(Np/2)-1);
   rangeGrid = (0:zeroPad)*Ts*c/2;
   
   %fig = figure('visible','off');
   fig = figure();
   set(fig,'Position',[0,30,1000,1000]);
   subplot(2,2,1)
   imagesc(rangeGrid,1:64,abs(matcharrayTx1));
   colorbar
   xlabel('Range (meters)'); ylabel('Pulse No');
   title(append('Delay vs Pulse Map TX1 CPI ',int2str(CPI)));
   
   subplot(2,2,3)
   imagesc(rangeGrid,nugrid, abs(rangedopplerTx1));
   xlabel('Range (meters)'); ylabel('Normalized Frequency');
   title(append('Doppler Range Map TX1 CPI ',int2str(CPI)));
   colorbar
   
   subplot(2,2,2)
   imagesc(rangeGrid,1:64,abs(matcharrayTx2));
   xlabel('Range (meters)'); ylabel('Pulse No');
   title(append('Delay vs Pulse Map TX2 CPI ',int2str(CPI)));
   colorbar
   
   subplot(2,2,4)
   imagesc(rangeGrid,nugrid, abs(rangedopplerTx2));
   xlabel('Range (meters)'); ylabel('Normalized Frequency');
   title(append('Doppler Range Map TX2 CPI ',int2str(CPI)));
   colorbar
   
   exportgraphics(fig,append('RDMap_CPI_',int2str(CPI),'.pdf'))
       
end