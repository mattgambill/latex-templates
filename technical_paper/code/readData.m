clear
fid=fopen('measured_data.bin','rb');
%aa = fread(fid,'double');
M=256; % number of pulses in one CPI
N_cpi= 32; % number of CPIs in data

 
for ii=1:N_cpi
    phase_history1(:,:,ii) = zeros(301,M);
    phase_history2(:,:,ii) = zeros(301,M);
    
    for jj=1:M
      phase_history1(:,jj,ii)=fread(fid,301,'double');
      phase_history1(:,jj,ii)=phase_history1(:,jj)+1i*fread(fid,301,'double');

      phase_history2(:,jj,ii)=fread(fid,301,'double');
      phase_history2(:,jj,ii)=phase_history2(:,jj)+1i*fread(fid,301,'double');


  
    end
  
   
   % only Range bins >=131 have target information, earlier range bins are
   % due to antenna coupling, TX leaking into RX, delays in the radar.
   % row range-bins 0:170, columns :pulse_no 1:M
   
   
    
   % phase_history matrix: row range-bins 0:170, columns :pulse_no 1:M
   % Construct Range Doppler map from the two antennas and display.
   % Consider canceling stationary clutter (two/three pulse canceller, or 
   % subtracting average (across slow time) range-profile from each return 
   
     
end
fclose(fid);
for ii=1:N_cpi
    phase_history1a(:,:,ii)=phase_history1(131:end,:,ii)';
    phase_history2a(:,:,ii)=phase_history2(131:end,:,ii)';
end
%%
disp('Done Reading')
disp('Writing Plots')
nugrid=1/M*(-M/2:1:(M/2)-1);
taugrid = (1:256).*(1/(250e6)).*299792458/2;

for ii=1:N_cpi

    fig = figure('visible','off');
    set(fig,'Position',[0,30,1000,1000]);
    rangedopplerTx1 = fftshift( fft(phase_history1a(:,:,ii),[],1),1);
    rangedopplerTx2 = fftshift( fft(phase_history2a(:,:,ii),[],1),1);
    
    subplot(2,2,1)
    imagesc(taugrid,1:256,abs(phase_history1a(:,:,ii)))
    title(append('Tx1 Range-Pulse Map CPI ',int2str(ii)))
    xlabel('Range (meters)');
    ylabel('Pulse Number');
    
    subplot(2,2,3)
    imagesc(taugrid,nugrid,abs(rangedopplerTx1))
    title(append('Tx1 Range-Doppler Map CPI ',int2str(ii)))
    xlabel('Range (meters)');
    ylabel('Normalized Doppler');
    
    
    subplot(2,2,2)
    imagesc(taugrid,1:256,abs(phase_history2a(:,:,ii)))
    title(append('Tx2 Range-Pulse Map CPI ',int2str(ii)))
    xlabel('Range (meters)');
    ylabel('Pulse Number');
    
    subplot(2,2,4)
    imagesc(taugrid,nugrid,abs(rangedopplerTx2))
    title(append('Tx2 Range-Doppler Map CPI ',int2str(ii)))
    xlabel('Range (meters)');
    ylabel('Normalized Doppler');
    
    exportgraphics(fig,append('dataPlot_',int2str(ii),'.pdf'))
 end