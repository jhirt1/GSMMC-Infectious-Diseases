function [SARdays] = updateRecoveryDays(R,SARdays)

ind=find(R); % Find the recovered individuals
SARdays(ind) = SARdays(ind) + 1;