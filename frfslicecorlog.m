function [newfslice, newfrf1slice, newfrf2slice, lagfactorvec] = frfslicecorlog(frf1,frf2,f,scaling,prom,percentoverlap,plotting)
% Levi Manring, Duke University
% January 2022
% This function takes two complex FRFs and uses logarithmic frequency
% interpolation to shift frf2 with respect to frf1. It uses
% cross-correlation to find the most similarity between the FRFs and then
% performs the shift. In particular, this function slices the frfs to find
% the regions of greatest similarity
%
% Inputs:   frf1 - complex vector containing the baseline FRF
%           frf2 - complex vector containing the comparison FRF
%           f - a range vector containing the frequencies for frf1 and frf2
%           scaling - a string either 'log' or 'any' that chooses whether
%                   or not to scale the frfs before performing cross-correlation
%           prom - prominence used to select peaks/features of frf1
%           percentoverlap - a percentage (0-1) used to denote how much
%                   overlap to include after slicing the frfs
%           plotting - a boolean (true or false) choosing to plot or not
%
% Outputs:  newfslice - an updated range vector containing frequences for the
%                   new FRFs
%           newfrf1slice - an updated frf1 interpolated at the same points as
%                   newfrf2slice
%           newfrf2slice - a sliced and shifted frf2 in the log frequency space
%           lagfactorvec - a transformation applied to the original f at each slice 
%                   to get a shifted frf2: i.e., f_frf2_shifted = lagfactor*f_frf2
%
%%
[~, frf1locs] = findpeaks(abs(frf1),'MinPeakProminence',prom);

% slice the frfs
rangevec = zeros(length(frf1locs),2);
rangevec(1,1) = 1;
rangevec(end,2) = length(frf1);
plac = round(diff(frf1locs)/2 + frf1locs(1:end-1));
for k = 1:(length(frf1locs)-1)
    rangevec(k,2) = plac(k,1);
    rangevec(k+1,1) = plac(k,1) + 1;
end
% modify for overlap
rangediff = diff(rangevec')';
sampleoverlap = round(percentoverlap*rangediff);
rangevecmod = rangevec;
for k = 1:(length(frf1locs) - 1)
    rangevecmod(k,2) = rangevec(k,2) + sampleoverlap(k,1);
    rangevecmod(k+1,1) = rangevec(k+1,1) - sampleoverlap(k,1);
end

% calculate the correlation for each slice
totalnewf = [];
totalnewfrf2 = [];
totalnewfrf1 = [];
lagfactorvec = zeros(length(rangevecmod),1);
for k = 1:length(rangevecmod)
    fsub = f(rangevecmod(k,1):rangevecmod(k,2),1);
    frf1sub = frf1(rangevecmod(k,1):rangevecmod(k,2),1);
    frf2sub = frf2(rangevecmod(k,1):rangevecmod(k,2),1);
    [newfsub, newfrf1sub, newfrf2sub, lagfactor] = frfcorlog(frf1sub,frf2sub,fsub,scaling);
    lagfactorvec(k,1) = lagfactor;
    totalnewf = [totalnewf; newfsub];
    totalnewfrf1 = [totalnewfrf1; newfrf1sub];
    totalnewfrf2 = [totalnewfrf2; newfrf2sub];
end
% process overlap sections using averaging
[new_f, new_frfs] = consolidate(totalnewf, [totalnewfrf1, totalnewfrf2]);

% output vectors
newfslice = new_f;
newfrf1slice = new_frfs(:,1);
newfrf2slice = new_frfs(:,2);

% show process via plotting if user desires
if plotting == true
    figure;
    subplot(4,1,1);
    loglog(f,abs(frf1),'-'); hold on;
    plot(f,abs(frf2),'-');
    xlabel('Frequency (Hz)'); ylabel('|H|');
    title('Original FRFs');
    [newf, newfrf1, newfrf2, ~] = frfcorlog(frf1,frf2,f,scaling);
    subplot(4,1,2);
    loglog(newf,abs(newfrf1),'-'); hold on;
    plot(newf,abs(newfrf2),'-');
    xlabel('Frequency(Hz)'); ylabel('|H|');
    title('Single-Shifted FRFs');
    subplot(4,1,3);
    loglog(totalnewf,abs(totalnewfrf1),'-'); hold on;
    plot(totalnewf,abs(totalnewfrf2),'-'); 
    xlabel('Frequency(Hz)'); ylabel('|H|');
    title('Slice-Shifted FRFs w/ Overlap');
    subplot(4,1,4);
    loglog(newfslice,abs(newfrf1slice),'-'); hold on;
    plot(newfslice,abs(newfrf2slice),'-');
    xlabel('Frequency (Hz)'); ylabel('|H|');
    title('Slice-Shifted FRFs w/o Overlap');

    figure;
    subplot(4,1,1);
    semilogx(f,unwrap(angle(frf1)),'-'); hold on;
    plot(f,unwrap(angle(frf2)),'-');
    xlabel('Frequency (Hz)'); ylabel('H phase');
    title('Original FRFs');
    [newf, newfrf1, newfrf2, ~] = frfcorlog(frf1,frf2,f,scaling);
    subplot(4,1,2);
    semilogx(newf,unwrap(angle(newfrf1)),'-'); hold on;
    plot(newf,unwrap(angle(newfrf2)),'-');
    xlabel('Frequency (Hz)'); ylabel('H phase');
    title('Single-Shifted FRFs');
    subplot(4,1,3);
    semilogx(totalnewf,unwrap(angle(totalnewfrf1)),'-'); hold on;
    plot(totalnewf,unwrap(angle(totalnewfrf2)),'-'); 
    xlabel('Frequency (Hz)'); ylabel('H phase');
    title('Slice-Shifted FRFs w/ Overlap');
    subplot(4,1,4);
    semilogx(newfslice,unwrap(angle(newfrf1slice)),'-'); hold on;
    plot(newfslice,unwrap(angle(newfrf2slice)),'-');
    xlabel('Frequency (Hz)'); ylabel('H phase');
    title('Slice-Shifted FRFs w/o Overlap');
end

end

function [xnew, ynew] = consolidate(x,y)
% find duplicate values of abscissa coordinate x, then average those y
% values to obtain a single value at each unique x point
[uniquex, iux, ~]  = unique(x);
uniquey = y(iux,:);
indextodupes = find(not(ismember(1:numel(x),iux)));
for k = 1:length(indextodupes)
    ydupes = y(find(x == x(indextodupes(k),1)),:);
    uniquey(find(uniquex == x(indextodupes(k),1)),:) = mean(ydupes);
end
xnew = uniquex;
ynew = uniquey;
end




