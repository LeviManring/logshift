function [newf, newfrf1, newfrf2, lagfactor] = frfcorlog(frf1,frf2,f,scaling)
% Levi Manring, Duke University
% January 2022
% This function takes two complex FRFs and uses logarithmic frequency
% interpolation to shift frf2 with respect to frf1. It uses
% cross-correlation to find the most similarity between the FRFs and then
% performs the shift.
%
% Inputs:   frf1 - complex vector containing the baseline FRF
%           frf2 - complex vector containing the comparison FRF
%           f - a range vector containing the frequencies for frf1 and frf2
%           scaling - a string either 'log' or 'any' that chooses whether
%                   or not to scale the frfs before performing cross-correlation
%
% Outputs:  newf - an updated range vector containing frequences for the
%                   new FRFs
%           newfrf1 - an updated frf1 interpolated at the same points as
%                   newfrf2
%           newfrf2 - a shifted frf2 in the log frequency space
%           lagfactor - a transformation applied to the original f to get
%                   a shifted frf2: i.e., f_frf2_shifted = lagfactor*f_frf2
%
%%
if scaling == "log"
    magfrf1 = log(abs(frf1));
    magfrf2 = log(abs(frf2));
else
    magfrf1 = abs(frf1);
    magfrf2 = abs(frf2);
end
% get the unwrapped angles for each FRF
anglefrf1 = unwrap(angle(frf1));
anglefrf2 = unwrap(angle(frf2));

[newlogf, newlogstuff] = log_interp(f,[magfrf1, magfrf2, anglefrf1, anglefrf2],10);
[~, nonlogstuff] = log_interp(f,[abs(frf1), abs(frf2)],10);
[c,lags] = xcorr(newlogstuff(:,1),newlogstuff(:,2));
[~, ind] = max(c);
maxlags = lags(ind);

logftoshift = maxlags*min(diff(newlogf));
lagfactor = exp(logftoshift);
newlogf2 = newlogf + logftoshift;
newlogf1 = newlogf;

newf1 = exp(newlogf1);
newf2 = exp(newlogf2);
realind1 = isfinite(newf1);
realind2 = isfinite(newf2);
newfrf1 = nonlogstuff(realind1,1).*exp(1i.*newlogstuff(realind1,3));
newfrf2 = nonlogstuff(realind2,2).*exp(1i.*newlogstuff(realind2,4));

newfrf1 = interp1(newf1,newfrf1,f);
newfrf2 = interp1(newf2,newfrf2,f);

nan1 = find(isnan(newfrf1));
nan2 = find(isnan(newfrf2));
allnans = union(nan1,nan2);

newfrf1(allnans) = [];
newf = f;
newf(allnans) = [];
newfrf2(allnans) = [];
end

function [newlogx, newlogy] = log_interp(x,y,ord)
% This function interpolates y using a logarithmic x scale
xnew = (min(x):(min(diff(x))/ord):max(x))';
ynew = interp1(x,y,xnew);
logx = log(xnew);
logx = logx(isfinite(log(xnew)));
ynew = ynew(isfinite(log(xnew)),:);
mindx = min(diff(logx));
newlogx = (min(logx):mindx:max(logx))';
newlogy = interp1(logx,ynew,newlogx);
end
