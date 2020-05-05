function [] = Make_All(wave,Fc)
%% Make_All:
%   INPUTS:
%       wave    :   spatial wavelengths of data
%       Fc      :   head data low-pass filter cutoff frequency [Hz]
%
%   OUTPUTS:
%       -
%

for w = 1:length(wave)
   MakeData_Ramp_HeadFree_Sacd(wave(w), Fc)    
end

disp('----------------------------')
disp('ALL DONE')
end