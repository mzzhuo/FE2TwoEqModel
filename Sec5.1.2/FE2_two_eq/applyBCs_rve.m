if strcmp(bc,'linear') 
    applyLinearBCs;
elseif strcmp(bc,'periodic') 
    applyPeriodicBCs;
else
    error('NO boundary condition type specified!');
end