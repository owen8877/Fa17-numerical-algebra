function rou = rayleighQuotient(U, Ah, Mh)
    rou = (U'*Ah*U) / (U'*Mh*U);
end

