function [H, R, ImmunityTime] = InfectedToRecovered(H, R, Hdays, rcd, ImmunityTime)

ind = find(H);
newrec = Hdays(ind) >= rcd;

% Update health statuses
idx = ind(newrec);
R(idx) = 1;
H(idx) = 0;

if ~isempty(idx)
    p = 0.7;
    mu1 = 20;   % Fast waning group average immunity
    mu2 = 200;  % Slow waning group average immunity 
    
    u = rand(length(idx), 1);
    fast = u < p;
    
    % Use explicit sizing [num_elements, 1]
    if any(fast)
        ImmunityTime(idx(fast))  = exprnd(mu1, [sum(fast), 1]);
    end
    if any(~fast)
        ImmunityTime(idx(~fast)) = exprnd(mu2, [sum(~fast), 1]);
    end
end
end
