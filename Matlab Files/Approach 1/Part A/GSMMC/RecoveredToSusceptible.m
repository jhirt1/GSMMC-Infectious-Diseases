function [R,S,SARdays] = RecoveredToSusceptible(R,S,SARdays,ImmunityTime)

inr = find(R);

rec = SARdays(inr) >= ImmunityTime(inr);

idx = inr(rec);

S(idx) = 1;
R(idx) = 0;
SARdays(idx) = 0;

end