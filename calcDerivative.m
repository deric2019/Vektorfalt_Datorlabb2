function dT = calcDerivative(next,previous,timestep)
    dT = abs((next-previous)/timestep);
end
