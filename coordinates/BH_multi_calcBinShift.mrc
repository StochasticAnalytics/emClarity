function [ binShift, binShiftTemplateSearch ] = BH_multi_calcBinShift(coords, samplingRate)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

tmpShift = (coords-fix(coords./samplingRate).*samplingRate);
binShiftTemplateSearch = [tmpShift(1), tmpShift(3) + tmpShift(2),tmpShift(4)];
binShift = binShiftTemplateSearch ./ samplingRate;

end

