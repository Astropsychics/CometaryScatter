
clear all
close all

for i={'8P','Encke','ISON','IZ','LS4','MH','PanSTARRS'}

comet_name = i{1};

fg = fopen(strcat('../../Results/',comet_name,'/Chandra_gas_spectrum_',comet_name,'.dat')); 
dg = fscanf(fg,'%f %f',[2,inf]); 
dg = dg';
fd = fopen(strcat('../../Results/',comet_name,'/Chandra_dust_spectrum_',comet_name,'.dat'));
dd = fscanf(fd,'%f %f',[2,inf]);
dd = dd'; 
fO = fopen(strcat('../Observations/',comet_name,'_intensity_Chandra.dat'));
dO = fscanf(fO,'%f %f',[2,inf]);
dO = dO';
fclose('all');

Eg = dg(:,1); 
Ig = dg(:,2); 
Ed = dd(:,1);
Id = dd(:,2);
EO = dO(:,1);
IO = dO(:,2); 

figure
semilogy(Eg,Ig,'g', Ed,Id,'r', EO,IO,'ko', 'LineWidth',2.5)
set(gca,'FontSize',16)
title(strcat(comet_name,' Scattering Comparison'))
xlabel('Energy [keV]')
ylabel('Counts sec^{-1} keV^{-1}');
legend('Model','Obs','Location','NorthEast');
print('-depsc2',strcat('../Chandra_gasdust_comparison_',comet_name,'.eps'));

clear all
close all
 
end
