function isiDisCheck(T, dt, pPreSpike, seed)
rng(seed);
preSpike = binornd(1, pPreSpike, 1, round(T/dt));
isi = diff(find(preSpike>0)*dt);

histogram(isi/dt);

end