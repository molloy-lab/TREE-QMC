### R code from vignette source 'zipfr-tutorial.Rnw'

###################################################
### code chunk number 1: zipfr-tutorial.Rnw:48-49
###################################################
options(useFancyQuotes=FALSE)


###################################################
### code chunk number 2: zipfr-tutorial.Rnw:80-81 (eval = FALSE)
###################################################
## install.packages("zipfR")


###################################################
### code chunk number 3: zipfr-tutorial.Rnw:93-94
###################################################
library(zipfR)


###################################################
### code chunk number 4: zipfr-tutorial.Rnw:98-99 (eval = FALSE)
###################################################
## ?zipfR


###################################################
### code chunk number 5: zipfr-tutorial.Rnw:146-147
###################################################
ItaRi.spc


###################################################
### code chunk number 6: zipfr-tutorial.Rnw:171-172
###################################################
summary(ItaRi.spc)


###################################################
### code chunk number 7: zipfr-tutorial.Rnw:180-182
###################################################
N(ItaRi.spc)
V(ItaRi.spc)


###################################################
### code chunk number 8: zipfr-tutorial.Rnw:193-195
###################################################
Vm(ItaRi.spc, 1)
Vm(ItaRi.spc, 1:5)


###################################################
### code chunk number 9: zipfr-tutorial.Rnw:203-204
###################################################
Vm(ItaRi.spc, 1) / N(ItaRi.spc)


###################################################
### code chunk number 10: zipfr-tutorial.Rnw:211-213 (eval = FALSE)
###################################################
## plot(ItaRi.spc)
## plot(ItaRi.spc, log="x")


###################################################
### code chunk number 11: zipfr-tutorial.Rnw:218-219
###################################################
plot(ItaRi.spc)


###################################################
### code chunk number 12: zipfr-tutorial.Rnw:222-223
###################################################
plot(ItaRi.spc, log="x")


###################################################
### code chunk number 13: zipfr-tutorial.Rnw:239-240
###################################################
with(ItaRi.spc, plot(m, Vm, main="Frequency Spectrum"))


###################################################
### code chunk number 14: zipfr-tutorial.Rnw:309-310
###################################################
head(ItaRi.emp.vgc)


###################################################
### code chunk number 15: zipfr-tutorial.Rnw:324-325
###################################################
ItaRi.emp.vgc


###################################################
### code chunk number 16: zipfr-tutorial.Rnw:329-330
###################################################
summary(ItaRi.emp.vgc)


###################################################
### code chunk number 17: zipfr-tutorial.Rnw:339-340
###################################################
plot(ItaRi.emp.vgc, add.m=1)


###################################################
### code chunk number 18: zipfr-tutorial.Rnw:369-371
###################################################
ItaRi.bin.vgc <- vgc.interp(ItaRi.spc, N(ItaRi.emp.vgc), m.max=1)
head(ItaRi.bin.vgc)


###################################################
### code chunk number 19: zipfr-tutorial.Rnw:395-397
###################################################
plot(ItaRi.emp.vgc, ItaRi.bin.vgc, 
     legend=c("observed", "interpolated"))


###################################################
### code chunk number 20: zipfr-tutorial.Rnw:428-430
###################################################
V(ItaRi.emp.vgc)[N(ItaRi.emp.vgc) == 50000]
V(ItaRi.spc)


###################################################
### code chunk number 21: zipfr-tutorial.Rnw:478-479
###################################################
ItaRi.fzm <- lnre("fzm", ItaRi.spc, exact=FALSE)


###################################################
### code chunk number 22: zipfr-tutorial.Rnw:494-495
###################################################
summary(ItaRi.fzm)


###################################################
### code chunk number 23: zipfr-tutorial.Rnw:513-514
###################################################
ItaRi.fzm.spc <- lnre.spc(ItaRi.fzm, N(ItaRi.fzm))


###################################################
### code chunk number 24: zipfr-tutorial.Rnw:519-520
###################################################
plot(ItaRi.spc, ItaRi.fzm.spc, legend=c("observed", "fZM"))


###################################################
### code chunk number 25: zipfr-tutorial.Rnw:531-532
###################################################
ItaRi.fzm.vgc <- lnre.vgc(ItaRi.fzm, (1:100) * 28e3)


###################################################
### code chunk number 26: zipfr-tutorial.Rnw:549-551
###################################################
plot(ItaRi.emp.vgc, ItaRi.fzm.vgc, N0=N(ItaRi.fzm), 
     legend=c("observed", "fZM"))


###################################################
### code chunk number 27: zipfr-tutorial.Rnw:612-614
###################################################
set.seed(42)
ItaRi.sub.spc <- sample.spc(ItaRi.spc, N=700000)


###################################################
### code chunk number 28: zipfr-tutorial.Rnw:621-623
###################################################
ItaRi.sub.fzm <- lnre("fzm", ItaRi.sub.spc, exact=FALSE)
ItaRi.sub.fzm


###################################################
### code chunk number 29: zipfr-tutorial.Rnw:633-634
###################################################
ItaRi.sub.fzm.vgc <- lnre.vgc(ItaRi.sub.fzm, N=N(ItaRi.emp.vgc))


###################################################
### code chunk number 30: zipfr-tutorial.Rnw:644-646
###################################################
  plot(ItaRi.bin.vgc, ItaRi.sub.fzm.vgc, N0=N(ItaRi.sub.fzm),
       legend=c("interpolated", "fZM"))


###################################################
### code chunk number 31: zipfr-tutorial.Rnw:685-687
###################################################
V(ItaUltra.spc)
V(ItaRi.spc)


###################################################
### code chunk number 32: zipfr-tutorial.Rnw:693-695
###################################################
N(ItaUltra.spc)
N(ItaRi.spc)


###################################################
### code chunk number 33: zipfr-tutorial.Rnw:702-704
###################################################
ItaUltra.fzm <- lnre("fzm", ItaUltra.spc, exact=FALSE)
ItaUltra.ext.vgc <- lnre.vgc(ItaUltra.fzm, N(ItaRi.emp.vgc))


###################################################
### code chunk number 34: zipfr-tutorial.Rnw:709-710
###################################################
plot(ItaUltra.ext.vgc, ItaRi.bin.vgc, legend=c("ultra-", "ri-"))


###################################################
### code chunk number 35: zipfr-tutorial.Rnw:724-725 (eval = FALSE)
###################################################
## data(package="zipfR")


###################################################
### code chunk number 36: zipfr-tutorial.Rnw:729-730 (eval = FALSE)
###################################################
## ?ItaRi.spc


###################################################
### code chunk number 37: zipfr-tutorial.Rnw:787-788
###################################################
V(Brown100k.spc)


###################################################
### code chunk number 38: zipfr-tutorial.Rnw:797-799
###################################################
Vseen <- V(Brown100k.spc) - Vm(Brown100k.spc, 1)
Vseen


###################################################
### code chunk number 39: zipfr-tutorial.Rnw:805-806
###################################################
Vseen / V(Brown100k.spc)


###################################################
### code chunk number 40: zipfr-tutorial.Rnw:812-813
###################################################
Vm(Brown100k.spc, 1) / N(Brown100k.spc)


###################################################
### code chunk number 41: zipfr-tutorial.Rnw:843-845
###################################################
Brown100k.zm <- lnre("zm", Brown100k.spc)
Brown100k.zm


###################################################
### code chunk number 42: zipfr-tutorial.Rnw:862-863
###################################################
EV(Brown100k.zm, c(1e6, 10e6, 100e6))


###################################################
### code chunk number 43: zipfr-tutorial.Rnw:870-871
###################################################
Vseen / EV(Brown100k.zm, c(1e6, 10e6, 100e6))


###################################################
### code chunk number 44: zipfr-tutorial.Rnw:876-877
###################################################
1 - (Vseen / EV(Brown100k.zm, c(1e6, 10e6, 100e6)))


###################################################
### code chunk number 45: zipfr-tutorial.Rnw:885-887
###################################################
N(Brown.spc)
V(Brown.spc)


###################################################
### code chunk number 46: zipfr-tutorial.Rnw:891-892
###################################################
EV(Brown100k.zm, N(Brown.spc)) 


###################################################
### code chunk number 47: zipfr-tutorial.Rnw:899-901
###################################################
1 - (Vseen / V(Brown.spc))
1 - (Vseen / EV(Brown100k.zm, N(Brown.spc)))


###################################################
### code chunk number 48: zipfr-tutorial.Rnw:927-928
###################################################
Brown.zm.spc <- lnre.spc(Brown100k.zm, N(Brown.spc))


###################################################
### code chunk number 49: zipfr-tutorial.Rnw:939-940
###################################################
EV(Brown100k.zm, N(Brown.spc)) - Vseen 


###################################################
### code chunk number 50: zipfr-tutorial.Rnw:946-948
###################################################
sum(Vm(Brown.zm.spc, 1))
sum(Vm(Brown.zm.spc, 1:2))


###################################################
### code chunk number 51: zipfr-tutorial.Rnw:951-952
###################################################
sum(Vm(Brown.zm.spc, 1:6))


###################################################
### code chunk number 52: zipfr-tutorial.Rnw:961-963
###################################################
Noov.zm <- sum(Vm(Brown.zm.spc, 1:6) * (1:6))
Noov.zm


###################################################
### code chunk number 53: zipfr-tutorial.Rnw:967-968
###################################################
Noov.zm / N(Brown.spc)


###################################################
### code chunk number 54: zipfr-tutorial.Rnw:977-978
###################################################
V(Brown.spc) - Vseen


###################################################
### code chunk number 55: zipfr-tutorial.Rnw:983-984
###################################################
sum(Vm(Brown.spc, 1:13))


###################################################
### code chunk number 56: zipfr-tutorial.Rnw:989-992
###################################################
Noov.emp <- sum(Vm(Brown.spc, 1:13) * (1:13))
Noov.emp
Noov.emp / N(Brown.spc)


###################################################
### code chunk number 57: zipfr-tutorial.Rnw:1014-1015
###################################################
Brown10M.zm.spc <- lnre.spc(Brown100k.zm, 10e6)


###################################################
### code chunk number 58: zipfr-tutorial.Rnw:1026-1027
###################################################
sum(Vm(Brown10M.zm.spc, 1:18) * (1:18))


###################################################
### code chunk number 59: zipfr-tutorial.Rnw:1032-1033
###################################################
sum(Vm(Brown10M.zm.spc, 1:18))


###################################################
### code chunk number 60: zipfr-tutorial.Rnw:1039-1040
###################################################
 EV(Brown100k.zm, 10e6) - sum(Vm(Brown10M.zm.spc, 1:18))


