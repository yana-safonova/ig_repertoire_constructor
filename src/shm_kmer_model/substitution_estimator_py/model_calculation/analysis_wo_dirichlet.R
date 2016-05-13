before.flu = 3

# Trivial vs k-neighbour

trivial.minus.kn <- function(x, y) {
  1 - y$N / x$N
}

hist(mapply(trivial.minus.kn, GMC.IGK.trivial[1:before.flu], GMC.IGK.k_neighbour[1:before.flu]), freq=F)
hist(mapply(trivial.minus.kn, GMC.IGL.trivial[1:before.flu], GMC.IGL.k_neighbour[1:before.flu]), freq=F)
hist(mapply(trivial.minus.kn, GMC.IGH.trivial[1:before.flu], GMC.IGH.k_neighbour[1:before.flu]), freq=F)
hist(mapply(trivial.minus.kn, GMC.IG.trivial[1:before.flu], GMC.IG.k_neighbour[1:before.flu]), freq=F)

hist(mapply(trivial.minus.kn, GMC.IGK.trivial[-(1:before.flu)], GMC.IGK.k_neighbour[-(1:before.flu)]), freq=F)
hist(mapply(trivial.minus.kn, GMC.IGL.trivial[-(1:before.flu)], GMC.IGL.k_neighbour[-(1:before.flu)]), freq=F)
hist(mapply(trivial.minus.kn, GMC.IGH.trivial[-(1:before.flu)], GMC.IGH.k_neighbour[-(1:before.flu)]), freq=F)
hist(mapply(trivial.minus.kn, GMC.IG.trivial[-(1:before.flu)], GMC.IG.k_neighbour[-(1:before.flu)]), freq=F)



hist(mapply(trivial.minus.kn, IDO.IGK.trivial[1:before.flu], IDO.IGK.k_neighbour[1:before.flu]), freq=F)
hist(mapply(trivial.minus.kn, IDO.IGL.trivial[1:before.flu], IDO.IGL.k_neighbour[1:before.flu]), freq=F)
hist(mapply(trivial.minus.kn, IDO.IGH.trivial[1:before.flu], IDO.IGH.k_neighbour[1:before.flu]), freq=F)
hist(mapply(trivial.minus.kn, IDO.IG.trivial[1:before.flu],  IDO.IG.k_neighbour[1:before.flu]), freq=F)

hist(mapply(trivial.minus.kn, IDO.IGK.trivial[-(1:before.flu)], IDO.IGK.k_neighbour[-(1:before.flu)]), freq=F)
hist(mapply(trivial.minus.kn, IDO.IGL.trivial[-(1:before.flu)], IDO.IGL.k_neighbour[-(1:before.flu)]), freq=F)
hist(mapply(trivial.minus.kn, IDO.IGH.trivial[-(1:before.flu)], IDO.IGH.k_neighbour[-(1:before.flu)]), freq=F)
hist(mapply(trivial.minus.kn, IDO.IG.trivial[-(1:before.flu)],  IDO.IG.k_neighbour[-(1:before.flu)]), freq=F)



hist(mapply(trivial.minus.kn, FV.IGK.trivial[1:before.flu], FV.IGK.k_neighbour[1:before.flu]), freq=F)
hist(mapply(trivial.minus.kn, FV.IGL.trivial[1:before.flu], FV.IGL.k_neighbour[1:before.flu]), freq=F)
hist(mapply(trivial.minus.kn, FV.IGH.trivial[1:before.flu], FV.IGH.k_neighbour[1:before.flu]), freq=F)
hist(mapply(trivial.minus.kn, FV.IG.trivial[1:before.flu],  FV.IG.k_neighbour[1:before.flu]), freq=F)

hist(mapply(trivial.minus.kn, FV.IGK.trivial[-(1:before.flu)], FV.IGK.k_neighbour[-(1:before.flu)]), freq=F)
hist(mapply(trivial.minus.kn, FV.IGL.trivial[-(1:before.flu)], FV.IGL.k_neighbour[-(1:before.flu)]), freq=F)
hist(mapply(trivial.minus.kn, FV.IGH.trivial[-(1:before.flu)], FV.IGH.k_neighbour[-(1:before.flu)]), freq=F)
hist(mapply(trivial.minus.kn, FV.IG.trivial[-(1:before.flu)],  FV.IG.k_neighbour[-(1:before.flu)]), freq=F)



median(mapply(trivial.minus.kn, GMC.IGK.trivial, GMC.IGK.k_neighbour), na.rm=TRUE)
median(mapply(trivial.minus.kn, GMC.IGL.trivial, GMC.IGL.k_neighbour), na.rm=TRUE)
median(mapply(trivial.minus.kn, GMC.IGH.trivial, GMC.IGH.k_neighbour), na.rm=TRUE)
median(mapply(trivial.minus.kn, GMC.IG.trivial, GMC.IG.k_neighbour), na.rm=TRUE)



median(mapply(trivial.minus.kn, IDO.IGK.trivial, IDO.IGK.k_neighbour), na.rm=TRUE)
median(mapply(trivial.minus.kn, IDO.IGL.trivial, IDO.IGL.k_neighbour), na.rm=TRUE)
median(mapply(trivial.minus.kn, IDO.IGH.trivial, IDO.IGH.k_neighbour), na.rm=TRUE)
median(mapply(trivial.minus.kn, IDO.IG.trivial, IDO.IG.k_neighbour), na.rm=TRUE)



median(mapply(trivial.minus.kn, FV.IGK.trivial, FV.IGK.k_neighbour), na.rm=TRUE)
median(mapply(trivial.minus.kn, FV.IGL.trivial, FV.IGL.k_neighbour), na.rm=TRUE)
median(mapply(trivial.minus.kn, FV.IGH.trivial, FV.IGH.k_neighbour), na.rm=TRUE)
median(mapply(trivial.minus.kn, FV.IG.trivial, FV.IG.k_neighbour), na.rm=TRUE)



mean(mapply(trivial.minus.kn, GMC.IGK.trivial, GMC.IGK.k_neighbour), na.rm=TRUE)
mean(mapply(trivial.minus.kn, GMC.IGL.trivial, GMC.IGL.k_neighbour), na.rm=TRUE)
mean(mapply(trivial.minus.kn, GMC.IGH.trivial, GMC.IGH.k_neighbour), na.rm=TRUE)
mean(mapply(trivial.minus.kn, GMC.IG.trivial, GMC.IG.k_neighbour), na.rm=TRUE)



mean(mapply(trivial.minus.kn, IDO.IGK.trivial, IDO.IGK.k_neighbour), na.rm=TRUE)
mean(mapply(trivial.minus.kn, IDO.IGL.trivial, IDO.IGL.k_neighbour), na.rm=TRUE)
mean(mapply(trivial.minus.kn, IDO.IGH.trivial, IDO.IGH.k_neighbour), na.rm=TRUE)
mean(mapply(trivial.minus.kn, IDO.IG.trivial, IDO.IG.k_neighbour), na.rm=TRUE)



mean(mapply(trivial.minus.kn, FV.IGK.trivial, FV.IGK.k_neighbour), na.rm=TRUE)
mean(mapply(trivial.minus.kn, FV.IGL.trivial, FV.IGL.k_neighbour), na.rm=TRUE)
mean(mapply(trivial.minus.kn, FV.IGH.trivial, FV.IGH.k_neighbour), na.rm=TRUE)
mean(mapply(trivial.minus.kn, FV.IG.trivial, FV.IG.k_neighbour), na.rm=TRUE)



cor.mm.method <- function(df1, df2) {
  diag(cor(df1, df2, use="na.or.complete"))
}



boxplot(t(mapply(cor.mm.method, GMC.IGK.trivial, GMC.IGK.k_neighbour)))
boxplot(t(mapply(cor.mm.method, GMC.IGL.trivial, GMC.IGL.k_neighbour)))
boxplot(t(mapply(cor.mm.method, GMC.IGH.trivial, GMC.IGH.k_neighbour)))
boxplot(t(mapply(cor.mm.method, GMC.IG.trivial, GMC.IG.k_neighbour)))


GMC 8 --- бедный датасет. После сборки баркодов его объем упал.
Этот датасет обладает минимальной корреляцией trivial и k_neighbour.
Для него минимальный объем выравниваний на V-gen в IGL модели (10KB). 


boxplot(t(mapply(cor.mm.method, GMC.IGK.trivial[1:before.flu], GMC.IGK.k_neighbour[1:before.flu])))
boxplot(t(mapply(cor.mm.method, GMC.IGL.trivial[1:before.flu], GMC.IGL.k_neighbour[1:before.flu])))
boxplot(t(mapply(cor.mm.method, GMC.IGH.trivial[1:before.flu], GMC.IGH.k_neighbour[1:before.flu])))
boxplot(t(mapply(cor.mm.method, GMC.IG.trivial[1:before.flu], GMC.IG.k_neighbour[1:before.flu])))



boxplot(t(mapply(cor.mm.method, IDO.IGK.trivial, IDO.IGK.k_neighbour)))
boxplot(t(mapply(cor.mm.method, IDO.IGL.trivial, IDO.IGL.k_neighbour)))
boxplot(t(mapply(cor.mm.method, IDO.IGH.trivial, IDO.IGH.k_neighbour)))
boxplot(t(mapply(cor.mm.method, IDO.IG.trivial, IDO.IG.k_neighbour)))



boxplot(t(mapply(cor.mm.method, FV.IGK.trivial, FV.IGK.k_neighbour)))
boxplot(t(mapply(cor.mm.method, FV.IGL.trivial, FV.IGL.k_neighbour)))
boxplot(t(mapply(cor.mm.method, FV.IGH.trivial, FV.IGH.k_neighbour)))
boxplot(t(mapply(cor.mm.method, FV.IG.trivial, FV.IG.k_neighbour)))


FV 5 -- очень странный датасет. Для него отсутствуют raw_reads. Возможно, что он 454. Датасет не пригоден для построения модели.

Таким образом, есть (исключая, FV 5 и GMC 8) консистентность before, after flu.
И на всех датасетах и типах цепей корреляция больше 85 процентов.


# IGH vs IGK


boxplot(t(mapply(cor.mm.method, GMC.IGH.k_neighbour[1:before.flu], GMC.IGK.k_neighbour[1:before.flu])))
boxplot(t(mapply(cor.mm.method, IDO.IGH.k_neighbour[1:before.flu], IDO.IGK.k_neighbour[1:before.flu])))
boxplot(t(mapply(cor.mm.method, FV.IGH.k_neighbour [1:before.flu], FV.IGK.k_neighbour [1:before.flu])))

boxplot(t(mapply(cor.mm.method, GMC.IGH.k_neighbour[-(1:before.flu)], GMC.IGK.k_neighbour[-(1:before.flu)])))
boxplot(t(mapply(cor.mm.method, IDO.IGH.k_neighbour[-(1:before.flu)], IDO.IGK.k_neighbour[-(1:before.flu)])))
boxplot(t(mapply(cor.mm.method, FV.IGH.k_neighbour [-(1:before.flu)], FV.IGK.k_neighbour [-(1:before.flu)])))

boxplot(t(mapply(cor.mm.method, GMC.IGH.k_neighbour, GMC.IGK.k_neighbour)))
boxplot(t(mapply(cor.mm.method, IDO.IGH.k_neighbour, IDO.IGK.k_neighbour)))
boxplot(t(mapply(cor.mm.method, FV.IGH.k_neighbour, FV.IGK.k_neighbour)))


# IGH vs IGL


boxplot(t(mapply(cor.mm.method, GMC.IGH.k_neighbour[1:before.flu], GMC.IGL.k_neighbour[1:before.flu])))
boxplot(t(mapply(cor.mm.method, IDO.IGH.k_neighbour[1:before.flu], IDO.IGL.k_neighbour[1:before.flu])))
boxplot(t(mapply(cor.mm.method, FV.IGH.k_neighbour [1:before.flu], FV.IGL.k_neighbour [1:before.flu])))

boxplot(t(mapply(cor.mm.method, GMC.IGH.k_neighbour[-(1:before.flu)], GMC.IGL.k_neighbour[-(1:before.flu)])))
boxplot(t(mapply(cor.mm.method, IDO.IGH.k_neighbour[-(1:before.flu)], IDO.IGL.k_neighbour[-(1:before.flu)])))
boxplot(t(mapply(cor.mm.method, FV.IGH.k_neighbour [-(1:before.flu)], FV.IGL.k_neighbour [-(1:before.flu)])))

boxplot(t(mapply(cor.mm.method, GMC.IGH.k_neighbour, GMC.IGL.k_neighbour)))
boxplot(t(mapply(cor.mm.method, IDO.IGH.k_neighbour, IDO.IGL.k_neighbour)))
boxplot(t(mapply(cor.mm.method, FV.IGH.k_neighbour, FV.IGL.k_neighbour)))


# IGK vs IGL


boxplot(t(mapply(cor.mm.method, GMC.IGK.k_neighbour[1:before.flu], GMC.IGL.k_neighbour[1:before.flu])))
boxplot(t(mapply(cor.mm.method, IDO.IGK.k_neighbour[1:before.flu], IDO.IGL.k_neighbour[1:before.flu])))
boxplot(t(mapply(cor.mm.method, FV.IGK.k_neighbour [1:before.flu], FV.IGL.k_neighbour [1:before.flu])))

boxplot(t(mapply(cor.mm.method, GMC.IGK.k_neighbour[-(1:before.flu)], GMC.IGL.k_neighbour[-(1:before.flu)])))
boxplot(t(mapply(cor.mm.method, IDO.IGK.k_neighbour[-(1:before.flu)], IDO.IGL.k_neighbour[-(1:before.flu)])))
boxplot(t(mapply(cor.mm.method, FV.IGK.k_neighbour [-(1:before.flu)], FV.IGL.k_neighbour [-(1:before.flu)])))

boxplot(t(mapply(cor.mm.method, GMC.IGK.k_neighbour, GMC.IGL.k_neighbour)))
boxplot(t(mapply(cor.mm.method, IDO.IGK.k_neighbour, IDO.IGL.k_neighbour)))
boxplot(t(mapply(cor.mm.method, FV.IGK.k_neighbour, FV.IGL.k_neighbour)))


Conclusion: bad correlations (~.5--.65).

# IG* vs IG* within individuals


t(outer(GMC.IGH.k_neighbour[1:before.flu], GMC.IGH.k_neighbour[1:before.flu], "+"))



# Individuals vs individuals


boxplot(t(mapply(cor.mm.method, GMC.IGH.k_neighbour, IDO.IGH.k_neighbour)))
boxplot(t(mapply(cor.mm.method, GMC.IGK.k_neighbour, IDO.IGK.k_neighbour)))
boxplot(t(mapply(cor.mm.method, GMC.IGL.k_neighbour, IDO.IGL.k_neighbour)))
boxplot(t(mapply(cor.mm.method, GMC.IG.k_neighbour, IDO.IG.k_neighbour)))

boxplot(t(mapply(cor.mm.method, GMC.IGH.k_neighbour[1:before.flu], IDO.IGH.k_neighbour[1:before.flu])))
boxplot(t(mapply(cor.mm.method, GMC.IGK.k_neighbour[1:before.flu], IDO.IGK.k_neighbour[1:before.flu])))
boxplot(t(mapply(cor.mm.method, GMC.IGL.k_neighbour[1:before.flu], IDO.IGL.k_neighbour[1:before.flu])))
boxplot(t(mapply(cor.mm.method, GMC.IG.k_neighbour [1:before.flu], IDO.IG.k_neighbour [1:before.flu])))

boxplot(t(mapply(cor.mm.method, GMC.IGH.k_neighbour[-(1:before.flu)], IDO.IGH.k_neighbour[-(1:before.flu)])))
boxplot(t(mapply(cor.mm.method, GMC.IGK.k_neighbour[-(1:before.flu)], IDO.IGK.k_neighbour[-(1:before.flu)])))
boxplot(t(mapply(cor.mm.method, GMC.IGL.k_neighbour[-(1:before.flu)], IDO.IGL.k_neighbour[-(1:before.flu)])))
boxplot(t(mapply(cor.mm.method, GMC.IG.k_neighbour [-(1:before.flu)], IDO.IG.k_neighbour [-(1:before.flu)])))



boxplot(t(mapply(cor.mm.method, GMC.IGH.k_neighbour, FV.IGH.k_neighbour)))
boxplot(t(mapply(cor.mm.method, GMC.IGK.k_neighbour, FV.IGK.k_neighbour)))
boxplot(t(mapply(cor.mm.method, GMC.IGL.k_neighbour, FV.IGL.k_neighbour)))
boxplot(t(mapply(cor.mm.method, GMC.IG.k_neighbour, FV.IG.k_neighbour)))



boxplot(t(mapply(cor.mm.method, IDO.IGH.k_neighbour, FV.IGH.k_neighbour)))
boxplot(t(mapply(cor.mm.method, IDO.IGK.k_neighbour, FV.IGK.k_neighbour)))
boxplot(t(mapply(cor.mm.method, IDO.IGL.k_neighbour, FV.IGL.k_neighbour)))
boxplot(t(mapply(cor.mm.method, IDO.IG.k_neighbour, FV.IG.k_neighbour)))