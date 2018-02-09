%% load B_scan
% load('B_scan.mat');
% flag = 0;
% if flag == 0
%     B_scan_data = flipud(B_scan_data);
% end
% flag = 1;
%% show the image
figure(1);
imagesc(B_scan_data);colormap jet;axis off;axis equal
figure(2);
one_spec = B_scan_data(:,1);
ctr_spc = 568.5;    % in nm
bandwidth = 107;    % in nm
spec_axis = linspace(ctr_spc - bandwidth/2, ctr_spc + bandwidth/2,length(one_spec));
plot(spec_axis,one_spec);xlabel('wavelength (nm)');ylabel('pixel intensity');axis square

%% smooth the curve
sm_stp = 1;
B_process = zeros(size(B_scan_data));
spec = mean(B_scan_data,2);
spec_smooth = smooth(spec,sm_stp);
for i = 1:size(B_scan_data,2)
    spec = double(B_scan_data(:,i));
    B_process(:,i) = spec - spec_smooth;
%     figure;
%     subplot(121);
%     plot(spec_smooth);axis square
%     subplot(122);
%     plot(spec - spec_smooth);axis square
end
figure(3);
imagesc(B_process);colormap jet;axis off;axis equal
figure(4);
one_spec = B_process(:,1);
plot(spec_axis,one_spec);xlabel('wavelength (nm)');ylabel('pixel intensity');axis square

%% get the corresponding frequency   
% % in Fourier domain zero padding to increase accuracy
% for i = 1:size(B_scan_data,2)
%     tmp = B_process(:,i);
%     tmp_f = zeros(2*length(tmp),1);
%     tmp_f(1:length(tmp)) = fftshift(fft(tmp));
%     tmp = ifft(ifftshift(tmp_f));
%     B2(:,i) = tmp;
% end
% spec_axis = linspace(ctr_spc - bandwidth/2, ctr_spc + bandwidth/2,size(B2,1));
frequency = 2*pi./spec_axis;
samp_num = 2*length(frequency);
min_f = min(frequency);
max_f = max(frequency);
freq_uni = linspace(min_f,max_f,samp_num);
% zero_num = round(min_f/(max_f-min_f)*samp_num);
% OCT = zeros(length(freq_uni),size(B_scan_data,2));
clear OCT
for i = 1:size(B_scan_data,2)
    spec_uni = interp1(frequency,B2(:,i),freq_uni,'spline');
%     for j = 1:samp_num
%         spec_zp(zero_num+j) = spec_uni(j);
%     end
    oct = fftshift(ifft(ifftshift(spec_uni)));
    OCT(:,i) = oct;
end
OCT_final = OCT(1:end/2,:);
figure(5);
imagesc(abs(OCT_final));axis off;axis equal
figure(6);
one_spec = OCT_final(:,1);
plot(abs(one_spec)*1e4);axis([0 2034 0 4e4]);