source('utils.r')

X = list(
SCRIPT = 'sim/app/mass.r',
NAME   = 'mass',
FLAGS  = '.debug=0 .k=RRo.rev.eR2',
NR     = 3000,   # num runs (total)
TPR    = 20,     # time per run (sec)
MPR    = 2,      # mem per run (GB)
NB     = 100,    # num batches
CPU    = 1)      # num cpu (parallel)

X$HMS = str(hms::hms(X$NR * X$TPR / X$NB / X$CPU * 2))
X$MEM = ceiling(X$CPU * X$MPR * 2)

status(2,list.str(X,join='\n'))
txt = load.txt('code','hpc','template',ext='.sh')
for (k in names(X)){
  txt = sub(str('<<',k,'>>'),X[[k]],txt) }
save.txt(txt,'code','hpc',X$NAME,ext='.sh')

