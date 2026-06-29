function [R,S,SARdays] = RecoveredToSusceptible_V(R,S,SARdays,ImmunityTime)
inr = find(R); % find the recovered individuals
rec = SARdays(inr) >= ImmunityTime(inr); % check how long they have stayed recovered
idx = inr(rec);
S(idx) = 1;
R(idx) = 0;
SARdays(idx) = 0;
end