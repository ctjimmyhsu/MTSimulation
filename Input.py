##################################################

#7.5 uM

tottime = 20*60 #minutes
interval = 5 #seconds

vg = 0.540477 #um/min
vs = 15 #um/min

bypassflag = 1 #bypass nucleation lag
nuctime = 200 #sec

numtubes = 10000

steps = 1.82185182429688
rate = 0.73924322088006733/60
loc = 0 #0.43750529414668771*60


lifetimes = simulate(tottime, interval, vg, vs, numtubes, steps, rate, loc, nuctime, bypassflag)