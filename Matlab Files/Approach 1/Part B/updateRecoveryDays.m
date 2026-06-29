function SARdays = updateRecoveryDays(R, SARdays)
    % Increment days for currently recovered individuals
    SARdays(R == 1) = SARdays(R == 1) + 1;
end