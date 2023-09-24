setwd('z:/2020/Raj_Paper/CleanData/blast_orfs/')

Depth<-read.table(file = 'z:/2020/Raj_Paper/CleanData/blast_orfs/BCTV2016_forR.fa', sep='\t',header = F)
Depth1<-read.table(file = 'z:/2020/Raj_Paper/CleanData/blast_orfs/BCTV2016_forR1.fa', sep='\t',header = F)
DepthMD<-read.table(file = 'z:/2020/Raj_Paper/CleanData/blast_orfs/BCTV_metadata.txt', sep=' ',header = F)


cuts <- function(x)
{
  n <- length(x) %/% 4
  map <- rep(c(rep(TRUE,4),FALSE), n)
  result <- rep(NA, n*5)
  result[map] <- x
  result
}


DepthMD

RN<-Depth[,1]
rownames(DepthMD)<-RN
DepthMD
rownames(Depth1)<-RN
Depth1<-Depth1[,-1]
Depth1<-as.matrix(Depth1)

Depth1[Depth1 == 0] <- 5
Depth1[Depth1 == 'M'] <- 5
Depth1[Depth1 == 'N'] <- 5
Depth1[Depth1 == 'R'] <- 5
Depth1[Depth1 == 'S'] <- 5
Depth1[Depth1 == 'Y'] <- 5

#dist(Depth1)
#table(Depth1)
Depth1<-as.matrix(Depth1)
pheatmap(Depth1,cluster_cols=F,fontsize_row = 1,cutree_rows = 1,clustering_method = 'ward.D2')
#hist(Depth2)
#class(Depth1)
Depth1[,c(2,3,4,5,6)]
Depth1<-data.matrix(Depth1)
#table(as.numeric(Depth1))


Depth2<-matrix(as.numeric(Depth1),ncol = 3175)
rownames(Depth2)<-RN
#pheatmap(Depth2[,seq(1,500)],cluster_cols=F,fontsize_row = 5,labels_row = F,cutree_rows = 1,clustering_method = 'ward.D2')
pheatmap(Depth2,cluster_cols=F,labels_row = F,fontsize_row = 5,cutree_rows = 1,clustering_method = 'ward.D2')
#heatmap(Depth2[,seq(1,500)])

color=colorRampPalette(c("green",'orange','lightblue','red','white')(5))

d <- dist(Depth2) # euclidean distances between the rows
fit <- cmdscale(d,eig=TRUE, k=2) # k is the number of dim
fit # view results

# plot solution
x <- fit$points[,1]
y <- fit$points[,2]

SI<-as.data.frame(unique(as.numeric(DepthMD$V2)))
SI
orderN<-order(SI)
color_colors<-rainbow(15)

plot(x, y, xlab="Coordinate 1", ylab="Coordinate 2",
     main="BCTV Virus Genome Sequences", type="p",pch=19,col=color_colors[as.numeric(DepthMD$V2)])
legend(x = "topleft",          # Position
       legend = unique(DepthMD$V2)[order(SI)],  # Legend texts
       pch = 19,           # Line types
       col = color_colors,           # Line colors
       lwd = 0)                 # Line width

row.names(Depth2)
text(x+1.1, y+0.1, labels = seq(1,length(row.names(Depth2))), cex=.7,col = 'black') 

as.vector(Depth2[,3])
hclust(Depth2)
CM<-cor(t(Depth2), method = "pearson", use = "complete.obs")
CM[c(2,4,5),c(2,4,5)]
pheatmap(CM,cluster_cols=F,fontsize_row = 5,cutree_rows = 1,clustering_method = 'ward.D2')


unique(DepthMD[,2])
CALRN<-DepthMD[DepthMD[,2]=='BCTV-CALogan',]
WORCucurbitaRN<-DepthMD[DepthMD[,2]=='BTCV-Wor',]
WORRN<-DepthMD[DepthMD[,2]=='BCTV-Wor',]
SEVRN<-DepthMD[DepthMD[,2]=='BCTV-Svr',]
CORN<-DepthMD[DepthMD[,2]=='BCTV-CO',]
PeYDRN<-DepthMD[DepthMD[,2]=='BCTV-PeYD',]
MLDRN<-DepthMD[DepthMD[,2]=='BCTV-Mld',]
SvrPepRN<-DepthMD[DepthMD[,2]=='BCTV-SvrPep',]
PeCTRN<-DepthMD[DepthMD[,2]=='BCTV-PeCT',]
SpCTRN<-DepthMD[DepthMD[,2]=='BCTV-SpCT',]
KIMRN<-DepthMD[DepthMD[,2]=='BCTV-Kim1',]
LH71TRN<-DepthMD[DepthMD[,2]=='BTCV-LH71',]
LH71CRN<-DepthMD[DepthMD[,2]=='BCTV-LH71',]

length(Depth2[1,])

FstALL<-data.frame()
Cor1<-data.frame()
AFV<-data.frame(1,2,3,4)
colnames(AFV)<-c('CA','WOR','SEV','CO')
for (i in seq(1,3175)){
  i<-i
  print(i)
  X<-as.vector(Depth2[rownames(CALRN),i])
  Xa<-(length(X[X==1])/length(X))^2
  Xt<-(length(X[X==2])/length(X))^2
  Xg<-(length(X[X==3])/length(X))^2
  Xc<-(length(X[X==4])/length(X))^2
  Xx<-(length(X[X==5])/length(X))^2
  SX<-1-sum(Xa,Xt,Xg,Xc,Xx)
  V1<-c(Xa,Xt,Xg,Xc,Xx)
  SUMX<-(SX)
  
  X1<-as.vector(Depth2[rownames(WORRN),i])
  Xa1<-(length(X1[X1==1])/length(X1))^2
  Xt1<-(length(X1[X1==2])/length(X1))^2
  Xg1<-(length(X1[X1==3])/length(X1))^2
  Xc1<-(length(X1[X1==4])/length(X1))^2
  Xx1<-(length(X1[X1==5])/length(X1))^2
  SX1<-1-sum(Xa1,Xt1,Xg1,Xc1,Xx1)
  SUMX1<-(SX1)
  FstSite<-((SX+SX1)-(SX1))/(SX+SX1)
  V2<-c(Xa1,Xt1,Xg1,Xc1,Xx1)
  C1<-cor(V1,V2)
  
  X2<-as.vector(Depth2[rownames(SEVRN),i])
  Xa2<-(length(X2[X2==1])/length(X2))^2
  Xt2<-(length(X2[X2==2])/length(X2))^2
  Xg2<-(length(X2[X2==3])/length(X2))^2
  Xc2<-(length(X2[X2==4])/length(X2))^2
  Xx2<-(length(X2[X2==5])/length(X2))^2
  SX2<-1-sum(Xa2,Xt2,Xg2,Xc2,Xx2)
  V3<-c(Xa2,Xt2,Xg2,Xc2,Xx2)
  
  X3<-as.vector(Depth2[rownames(CORN),i])
  Xa3<-(length(X3[X3==1])/length(X3))^2
  Xt3<-(length(X3[X3==2])/length(X3))^2
  Xg3<-(length(X3[X3==3])/length(X3))^2
  Xc3<-(length(X3[X3==4])/length(X3))^2
  Xx3<-(length(X3[X3==5])/length(X3))^2
  SX3<-1-sum(Xa3,Xt3,Xg3,Xc3,Xx3)
  V4<-c(Xa3,Xt3,Xg3,Xc3,Xx3)
  
  C1cw<-cor(V1,V2)
  C1cs<-cor(V1,V3)
  C1cco<-cor(V1,V4)
  
  C1wc<-cor(V2,V1)
  C1ws<-cor(V2,V3)
  C1wco<-cor(V2,V4)
  
  
  C1sc<-cor(V3,V1)
  C1sw<-cor(V3,V2)
  C1sco<-cor(V3,V4)
  
  C1coc<-cor(V4,V1)
  C1cow<-cor(V4,V2)
  C1cos<-cor(V4,V3)
  
  C1cor<-c(C1cw,C1cs,C1cco,C1wc,C1ws,C1wco,C1sc,C1sw,C1sco,C1coc,C1cow,C1cos)
  AFVc<-SX
  AFVw<-SX1
  AFVs<-SX2
  AFVco<-SX3
  AFV1<-c(SX,SX1,SX2,SX3)
  AFV<-rbind(AFV,AFV1)
  Cor1<-rbind(Cor1,C1cor)
}
AFV<-AFV[-1,]
Cor1<-(1-abs(Cor1))

##Avergae divergence between strain###
mean(rowMeans(Cor1[,c(1,2,3)]))
mean(rowMeans(Cor1[,c(4,5,6)]))
mean(rowMeans(Cor1[,c(7,8,9)]))
mean(rowMeans(Cor1[,c(10,11,12)]))


#LOS<-loess.smooth(seq(1,3175),Cor1[,1],span = 1/50,evaluation = 1000)
AFV[AFV==0]<-0.00001
AFV[AFV<0]<-0.00001
colMeans(AFV)



par(mfrow=c(4,1))
par(oma=c(3,0.5,3,0.5),mar=c(1,3,3,0.5))
par(bty = 'n') 

plot(seq(1,3175),AFV[,1],col='red',type='l',cex=.3,ylim=c(-1,0.7))
par(new=T)


VecCor<-Cor1[,1]
n<-1
for (i in VecCor){
  cta<-(.3*i)
  color_transparent <- adjustcolor('dodgerblue', alpha.f = cta)
  xbottom <- n
  ybottom <- -0.30
  xtop <- n+1
  ytop<- -0.1
  xpol <- c(xbottom,xbottom,xtop,xtop)
  ypol <- c(ybottom,ytop,ytop,ybottom)
  polygon(x=cuts(t(xpol)), y=cuts(t(ypol)), col = color_transparent, border =color_transparent)
  n<-n+1
}

VecCor<-Cor1[,2]
n<-1
for (i in VecCor){
  cta<-(.3*i)
  color_transparent <- adjustcolor('orange', alpha.f = cta)
  xbottom <- n
  ybottom <- -0.50
  xtop <- n+1
  ytop<- -0.3
  xpol <- c(xbottom,xbottom,xtop,xtop)
  ypol <- c(ybottom,ytop,ytop,ybottom)
  polygon(x=cuts(t(xpol)), y=cuts(t(ypol)), col = color_transparent, border =color_transparent)
  n<-n+1
}

VecCor<-Cor1[,3]
n<-1
for (i in VecCor){
  cta<-(.3*i)
  color_transparent <- adjustcolor('green', alpha.f = cta)
  xbottom <- n
  ybottom <- -0.70
  xtop <- n+1
  ytop<- -0.5
  xpol <- c(xbottom,xbottom,xtop,xtop)
  ypol <- c(ybottom,ytop,ytop,ybottom)
  polygon(x=cuts(t(xpol)), y=cuts(t(ypol)), col = color_transparent, border =color_transparent)
  n<-n+1
}

plot(seq(1,3175),AFV[,2],col='dodgerblue',type='l',cex=.3,ylim=c(-1,0.7))
par(new=T)


VecCor<-Cor1[,4]
n<-1
for (i in VecCor){
  cta<-(.3*i)
  color_transparent <- adjustcolor('red', alpha.f = cta)
  xbottom <- n
  ybottom <- -0.30
  xtop <- n+1
  ytop<- -0.1
  xpol <- c(xbottom,xbottom,xtop,xtop)
  ypol <- c(ybottom,ytop,ytop,ybottom)
  polygon(x=cuts(t(xpol)), y=cuts(t(ypol)), col = color_transparent, border =color_transparent)
  n<-n+1
}

VecCor<-Cor1[,5]
n<-1
for (i in VecCor){
  cta<-(.3*i)
  color_transparent <- adjustcolor('orange', alpha.f = cta)
  xbottom <- n
  ybottom <- -0.50
  xtop <- n+1
  ytop<- -0.3
  xpol <- c(xbottom,xbottom,xtop,xtop)
  ypol <- c(ybottom,ytop,ytop,ybottom)
  polygon(x=cuts(t(xpol)), y=cuts(t(ypol)), col = color_transparent, border =color_transparent)
  n<-n+1
}

VecCor<-Cor1[,6]
n<-1
for (i in VecCor){
  cta<-(.3*i)
  color_transparent <- adjustcolor('green', alpha.f = cta)
  xbottom <- n
  ybottom <- -0.70
  xtop <- n+1
  ytop<- -0.5
  xpol <- c(xbottom,xbottom,xtop,xtop)
  ypol <- c(ybottom,ytop,ytop,ybottom)
  polygon(x=cuts(t(xpol)), y=cuts(t(ypol)), col = color_transparent, border =color_transparent)
  n<-n+1
}


plot(seq(1,3175),AFV[,3],col='orange',type='l',cex=.3,ylim=c(-1,0.7))
par(new=T)


VecCor<-Cor1[,7]
n<-1
for (i in VecCor){
  cta<-(.3*i)
  color_transparent <- adjustcolor('red', alpha.f = cta)
  xbottom <- n
  ybottom <- -0.30
  xtop <- n+1
  ytop<- -0.1
  xpol <- c(xbottom,xbottom,xtop,xtop)
  ypol <- c(ybottom,ytop,ytop,ybottom)
  polygon(x=cuts(t(xpol)), y=cuts(t(ypol)), col = color_transparent, border =color_transparent)
  n<-n+1
}

VecCor<-Cor1[,8]
n<-1
for (i in VecCor){
  cta<-(.3*i)
  color_transparent <- adjustcolor('dodgerblue', alpha.f = cta)
  xbottom <- n
  ybottom <- -0.50
  xtop <- n+1
  ytop<- -0.3
  xpol <- c(xbottom,xbottom,xtop,xtop)
  ypol <- c(ybottom,ytop,ytop,ybottom)
  polygon(x=cuts(t(xpol)), y=cuts(t(ypol)), col = color_transparent, border =color_transparent)
  n<-n+1
}

VecCor<-Cor1[,9]
n<-1
for (i in VecCor){
  cta<-(.3*i)
  color_transparent <- adjustcolor('green', alpha.f = cta)
  xbottom <- n
  ybottom <- -0.70
  xtop <- n+1
  ytop<- -0.5
  xpol <- c(xbottom,xbottom,xtop,xtop)
  ypol <- c(ybottom,ytop,ytop,ybottom)
  polygon(x=cuts(t(xpol)), y=cuts(t(ypol)), col = color_transparent, border =color_transparent)
  n<-n+1
}

plot(seq(1,3175),AFV[,4],col='green',type='l',cex=.3,ylim=c(-1,0.7))
par(new=T)


VecCor<-Cor1[,10]
n<-1
for (i in VecCor){
  cta<-(.3*i)
  color_transparent <- adjustcolor('red', alpha.f = cta)
  xbottom <- n
  ybottom <- -0.30
  xtop <- n+1
  ytop<- -0.1
  xpol <- c(xbottom,xbottom,xtop,xtop)
  ypol <- c(ybottom,ytop,ytop,ybottom)
  polygon(x=cuts(t(xpol)), y=cuts(t(ypol)), col = color_transparent, border =color_transparent)
  n<-n+1
}

VecCor<-Cor1[,11]
n<-1
for (i in VecCor){
  cta<-(.3*i)
  color_transparent <- adjustcolor('dodgerblue', alpha.f = cta)
  xbottom <- n
  ybottom <- -0.50
  xtop <- n+1
  ytop<- -0.3
  xpol <- c(xbottom,xbottom,xtop,xtop)
  ypol <- c(ybottom,ytop,ytop,ybottom)
  polygon(x=cuts(t(xpol)), y=cuts(t(ypol)), col = color_transparent, border =color_transparent)
  n<-n+1
}

VecCor<-Cor1[,12]
n<-1
for (i in VecCor){
  cta<-(.3*i)
  color_transparent <- adjustcolor('orange', alpha.f = cta)
  xbottom <- n
  ybottom <- -0.70
  xtop <- n+1
  ytop<- -0.5
  xpol <- c(xbottom,xbottom,xtop,xtop)
  ypol <- c(ybottom,ytop,ytop,ybottom)
  polygon(x=cuts(t(xpol)), y=cuts(t(ypol)), col = color_transparent, border =color_transparent)
  n<-n+1
}









####ORFs#####
{
  ORFs<-read.table('out.cds.blast')
  
  ORFF<-ORFs[ORFs[,1]=='giCAL',c(13,14)]
  
  
  
  lorf<-c()
  for (p in seq(1,length(ORFF[,1]))){
    ll<-(length(seq(ORFF[p,1],ORFF[p,2])))
    lorf<-c(lorf,ll)
    
  }
  lorf
  lorfo<-order(lorf,decreasing = F)
  
  Pdf<-c(rep(0,3175))
  Pfd1<-c(0,0,0)
  for (p in lorfo){
    print(ORFF[p,])
    P1<-c(rep(0,3175))
    P1[seq(ORFF[p,1],ORFF[p,2])]=1
    Pdf<-rbind(Pdf,P1)
    print((colSums(Pdf[,seq(ORFF[p,1],ORFF[p,2])])))
    MCS<-max(colSums(Pdf[,seq(ORFF[p,1],ORFF[p,2])]))
    ORFFF<-(cbind(ORFF[p,],MCS))
    Pfd1<-rbind(Pfd1,ORFFF)
  }
  
  Pfd1<-Pfd1[-1,]
  
  #plot(NULL,ylim=c(0,40),xlim = c(0,3175))
  par(oma=c(3,0.5,3,0.5),mar=c(1,3,3,0.5))
  par(bty = 'n') 
  plot(NULL,ylim=c(-2,40),xlim = c(0,3175))
  
  for (p in seq(1,length(Pfd1[,1]))){
    color_transparent <- adjustcolor('red', alpha.f = .5)
    xbottom <- Pfd1[p,1]
    ybottom <- Pfd1[p,3]+35
    xtop <- Pfd1[p,2]
    ytop<- Pfd1[p,3]+36
    xpol <- c(xbottom,xbottom,xtop,xtop)
    ypol <- c(ybottom,ytop,ytop,ybottom)
    polygon(x=cuts(t(xpol)), y=cuts(t(ypol)), col = color_transparent, border ='black')
  }
  
  
  ORFF<-ORFs[ORFs[,1]=='giWOR',c(13,14)]
  
  
  lorf<-c()
  for (p in seq(1,length(ORFF[,1]))){
    ll<-(length(seq(ORFF[p,1],ORFF[p,2])))
    lorf<-c(lorf,ll)
    
  }
  lorf
  lorfo<-order(lorf,decreasing = F)
  
  Pdf<-c(rep(0,3175))
  Pfd1<-c(0,0,0)
  for (p in lorfo){
    print(ORFF[p,])
    P1<-c(rep(0,3175))
    P1[seq(ORFF[p,1],ORFF[p,2])]=1
    Pdf<-rbind(Pdf,P1)
    print((colSums(Pdf[,seq(ORFF[p,1],ORFF[p,2])])))
    MCS<-max(colSums(Pdf[,seq(ORFF[p,1],ORFF[p,2])]))
    ORFFF<-(cbind(ORFF[p,],MCS))
    Pfd1<-rbind(Pfd1,ORFFF)
  }
  
  Pfd1<-Pfd1[-1,]
  
  #plot(NULL,ylim=c(0,40),xlim = c(0,3175))
  
  for (p in seq(1,length(Pfd1[,1]))){
    color_transparent <- adjustcolor('dodgerblue', alpha.f = .5)
    xbottom <- Pfd1[p,1]
    ybottom <- Pfd1[p,3]+27
    xtop <- Pfd1[p,2]
    ytop<- Pfd1[p,3]+28
    xpol <- c(xbottom,xbottom,xtop,xtop)
    ypol <- c(ybottom,ytop,ytop,ybottom)
    polygon(x=cuts(t(xpol)), y=cuts(t(ypol)), col = color_transparent, border ='black')
  }
  
  
  
  ORFF<-ORFs[ORFs[,1]=='giSEV',c(13,14)]
  
  lorf<-c()
  for (p in seq(1,length(ORFs[,1]))){
    ll<-length(seq(ORFs[p,13],ORFs[p,14]))
    lorf<-c(lorf,ll)
  }
  
  lorf<-c()
  for (p in seq(1,length(ORFF[,1]))){
    ll<-(length(seq(ORFF[p,1],ORFF[p,2])))
    lorf<-c(lorf,ll)
    
  }
  lorf
  lorfo<-order(lorf,decreasing = F)
  
  Pdf<-c(rep(0,3175))
  Pfd1<-c(0,0,0)
  for (p in lorfo){
    print(ORFF[p,])
    P1<-c(rep(0,3175))
    P1[seq(ORFF[p,1],ORFF[p,2])]=1
    Pdf<-rbind(Pdf,P1)
    print((colSums(Pdf[,seq(ORFF[p,1],ORFF[p,2])])))
    MCS<-max(colSums(Pdf[,seq(ORFF[p,1],ORFF[p,2])]))
    print(MCS)
    ORFFF<-(cbind(ORFF[p,],MCS))
    Pfd1<-rbind(Pfd1,ORFFF)
  }
  
  Pfd1<-Pfd1[-1,]
  
  
  for (p in seq(1,length(Pfd1[,1]))){
    color_transparent <- adjustcolor('orange', alpha.f = .5)
    xbottom <- Pfd1[p,1]
    ybottom <- Pfd1[p,3]+19
    xtop <- Pfd1[p,2]
    ytop<- Pfd1[p,3]+20
    xpol <- c(xbottom,xbottom,xtop,xtop)
    ypol <- c(ybottom,ytop,ytop,ybottom)
    polygon(x=cuts(t(xpol)), y=cuts(t(ypol)), col = color_transparent, border ='black')
  }
  
  
  
  
  
  ORFF<-ORFs[ORFs[,1]=='giCO',c(13,14)]
  
  
  lorf<-c()
  for (p in seq(1,length(ORFF[,1]))){
    ll<-(length(seq(ORFF[p,1],ORFF[p,2])))
    lorf<-c(lorf,ll)
    
  }
  lorf
  lorfo<-order(lorf,decreasing = F)
  
  Pdf<-c(rep(0,3175))
  Pfd1<-c(0,0,0)
  for (p in lorfo){
    print(ORFF[p,])
    P1<-c(rep(0,3175))
    P1[seq(ORFF[p,1],ORFF[p,2])]=1
    Pdf<-rbind(Pdf,P1)
    print((colSums(Pdf[,seq(ORFF[p,1],ORFF[p,2])])))
    MCS<-max(colSums(Pdf[,seq(ORFF[p,1],ORFF[p,2])]))
    ORFFF<-(cbind(ORFF[p,],MCS))
    Pfd1<-rbind(Pfd1,ORFFF)
  }
  
  Pfd1<-Pfd1[-1,]
  
  #plot(NULL,ylim=c(0,40),xlim = c(0,3175))
  
  for (p in seq(1,length(Pfd1[,1]))){
    color_transparent <- adjustcolor('green', alpha.f = .5)
    xbottom <- Pfd1[p,1]
    ybottom <- Pfd1[p,3]+11
    xtop <- Pfd1[p,2]
    ytop<- Pfd1[p,3]+12
    xpol <- c(xbottom,xbottom,xtop,xtop)
    ypol <- c(ybottom,ytop,ytop,ybottom)
    polygon(x=cuts(t(xpol)), y=cuts(t(ypol)), col = color_transparent, border ='black')
  }
  
  
  #plot(NULL,ylim=c(0,40),xlim = c(0,3175))
  sncRNA<-read.table('out.snc.blast')
  Pfd1<-sncRNA[,c(9,10)]
  for (p in seq(1,length(Pfd1[,1]))){
    color_transparent <- adjustcolor('black', alpha.f = .5)
    xbottom <- Pfd1[p,1]
    ybottom <- -2
    xtop <- Pfd1[p,2]
    ytop<- -1
    xpol <- c(xbottom,xbottom,xtop,xtop)
    ypol <- c(ybottom,ytop,ytop,ybottom)
    polygon(x=cuts(t(xpol)), y=cuts(t(ypol)), col = color_transparent, border ='black')
  }
}

#### SNCRNA strain specific
{
  sncRNA<-read.table('SNC_strainspecific.txt',sep = '\t')
  sncRNAT<-sncRNA
  sncRNAS<-sncRNA[sncRNA[,5]=='x',]
  sncRNACA<-sncRNA[sncRNA[,6]=='x',]
  sncRNATCO<-sncRNA[sncRNA[,7]=='x',]
  sncRNATWOR<-sncRNA[sncRNA[,8]=='x',]
  
  
  
  sncRNA<-read.table('out.snc.blast')
  sncRNA<-sncRNA[sncRNA[,1]%in%sncRNAS[,2]==T,]
  Pfd1<-sncRNA[,c(9,10)]
  for (p in seq(1,length(Pfd1[,1]))){
    color_transparent <- adjustcolor('black', alpha.f = .5)
    xbottom <- Pfd1[p,1]
    ybottom <- 19
    xtop <- Pfd1[p,2]
    ytop<- 20
    xpol <- c(xbottom,xbottom,xtop,xtop)
    ypol <- c(ybottom,ytop,ytop,ybottom)
    polygon(x=cuts(t(xpol)), y=cuts(t(ypol)), col = color_transparent, border ='black')
  }
  
  sncRNA<-read.table('out.snc.blast')
  sncRNA<-sncRNA[sncRNA[,1]%in%sncRNACA[,2]==T,]
  Pfd1<-sncRNA[,c(9,10)]
  for (p in seq(1,length(Pfd1[,1]))){
    color_transparent <- adjustcolor('black', alpha.f = .5)
    xbottom <- Pfd1[p,1]
    ybottom <- 35
    xtop <- Pfd1[p,2]
    ytop<- 36
    xpol <- c(xbottom,xbottom,xtop,xtop)
    ypol <- c(ybottom,ytop,ytop,ybottom)
    polygon(x=cuts(t(xpol)), y=cuts(t(ypol)), col = color_transparent, border ='black')
  }
  
  sncRNA<-read.table('out.snc.blast')
  sncRNA<-sncRNA[sncRNA[,1]%in%sncRNATWOR[,2]==T,]
  Pfd1<-sncRNA[,c(9,10)]
  for (p in seq(1,length(Pfd1[,1]))){
    color_transparent <- adjustcolor('black', alpha.f = .5)
    xbottom <- Pfd1[p,1]
    ybottom <- 27
    xtop <- Pfd1[p,2]
    ytop<- 28
    xpol <- c(xbottom,xbottom,xtop,xtop)
    ypol <- c(ybottom,ytop,ytop,ybottom)
    polygon(x=cuts(t(xpol)), y=cuts(t(ypol)), col = color_transparent, border ='black')
  }
  
  sncRNA<-read.table('out.snc.blast')
  sncRNA<-sncRNA[sncRNA[,1]%in%sncRNATCO[,2]==T,]
  Pfd1<-sncRNA[,c(9,10)]
  for (p in seq(1,length(Pfd1[,1]))){
    color_transparent <- adjustcolor('black', alpha.f = .5)
    xbottom <- Pfd1[p,1]
    ybottom <- 11
    xtop <- Pfd1[p,2]
    ytop<- 12
    xpol <- c(xbottom,xbottom,xtop,xtop)
    ypol <- c(ybottom,ytop,ytop,ybottom)
    polygon(x=cuts(t(xpol)), y=cuts(t(ypol)), col = color_transparent, border ='black')
  }
}

#BCTV genes
{
  sncRNA<-read.table('BCTVgenes_prot.blast.out')
  Gm<-c(min(sncRNA[sncRNA[,1]=='V1',c(9,10)]))
  GM<-c(max(sncRNA[sncRNA[,1]=='V1',c(9,10)]))
  xbottom <- Gm
  ybottom <- 5
  xtop <- GM
  ytop<- 6
  xpol <- c(xbottom,xbottom,xtop,xtop)
  ypol <- c(ybottom,ytop,ytop,ybottom)
  polygon(x=cuts(t(xpol)), y=cuts(t(ypol)), col = color_transparent, border ='black')
  
  Gm<-c(min(sncRNA[sncRNA[,1]=='V2',c(9,10)]))
  GM<-c(max(sncRNA[sncRNA[,1]=='V2',c(9,10)]))
  xbottom <- Gm
  ybottom <- 6
  xtop <- GM
  ytop<- 7
  xpol <- c(xbottom,xbottom,xtop,xtop)
  ypol <- c(ybottom,ytop,ytop,ybottom)
  polygon(x=cuts(t(xpol)), y=cuts(t(ypol)), col = color_transparent, border ='black')
  
  Gm<-c(min(sncRNA[sncRNA[,1]=='V3',c(9,10)]))
  GM<-c(max(sncRNA[sncRNA[,1]=='V3',c(9,10)]))
  xbottom <- Gm
  ybottom <- 5
  xtop <- GM
  ytop<- 6
  xpol <- c(xbottom,xbottom,xtop,xtop)
  ypol <- c(ybottom,ytop,ytop,ybottom)
  polygon(x=cuts(t(xpol)), y=cuts(t(ypol)), col = color_transparent, border ='black')
  
  Gm<-c(min(sncRNA[sncRNA[,1]=='C1',c(9,10)]))
  GM<-c(max(sncRNA[sncRNA[,1]=='C1',c(9,10)]))
  xbottom <- Gm
  ybottom <- 5
  xtop <- GM
  ytop<- 6
  xpol <- c(xbottom,xbottom,xtop,xtop)
  ypol <- c(ybottom,ytop,ytop,ybottom)
  polygon(x=cuts(t(xpol)), y=cuts(t(ypol)), col = color_transparent, border ='black')
  
  Gm<-c(min(sncRNA[sncRNA[,1]=='C2',c(9,10)]))
  GM<-c(max(sncRNA[sncRNA[,1]=='C2',c(9,10)]))
  xbottom <- Gm
  ybottom <- 6
  xtop <- GM
  ytop<- 7
  xpol <- c(xbottom,xbottom,xtop,xtop)
  ypol <- c(ybottom,ytop,ytop,ybottom)
  polygon(x=cuts(t(xpol)), y=cuts(t(ypol)), col = color_transparent, border ='black')
  
  sncRNA[sncRNA[,1]=='C3',c(9,10)]
  Gm<-c(min(sncRNA[sncRNA[,1]=='C3',c(9,10)]))
  GM<-c(max(sncRNA[sncRNA[,1]=='C3',c(9,10)]))
  xbottom <- Gm
  ybottom <- 5
  xtop <- GM
  ytop<- 6
  xpol <- c(xbottom,xbottom,xtop,xtop)
  ypol <- c(ybottom,ytop,ytop,ybottom)
  polygon(x=cuts(t(xpol)), y=cuts(t(ypol)), col = color_transparent, border ='black')
  
  sncRNA[sncRNA[,1]=='C4',c(9,10)]
  Gm<-c(min(sncRNA[sncRNA[,1]=='C4',c(9,10)]))
  GM<-c(max(sncRNA[sncRNA[,1]=='C4',c(9,10)]))
  xbottom <- Gm
  ybottom <- 6
  xtop <- GM
  ytop<- 7
  xpol <- c(xbottom,xbottom,xtop,xtop)
  ypol <- c(ybottom,ytop,ytop,ybottom)
  polygon(x=cuts(t(xpol)), y=cuts(t(ypol)), col = color_transparent, border ='black')
}

AFVrm<-rowMeans(AFV)
for (p in seq(1,length(AFVrm))){
  cta=AFVrm[p]
  color_transparent <- adjustcolor('white', alpha.f = cta)
  xbottom <- p
  ybottom <- 0
  xtop <- p
  ytop<- 1
  xpol <- c(xbottom,xbottom,xtop,xtop)
  ypol <- c(ybottom,ytop,ytop,ybottom)
  polygon(x=cuts(t(xpol)), y=cuts(t(ypol)), col = color_transparent, border = color_transparent)
}



DivCAL<-rowMeans(Cor1[,c(1,2,3)])
DivWOR<-rowMeans(Cor1[,c(4,5,6)])
DivSEV<-rowMeans(Cor1[,c(7,8,9)])
DivCO<-rowMeans(Cor1[,c(10,11,12)])


AFV
LOS1<-loess.smooth(seq(1,length(DivCAL)),DivCAL,degree = 1,
                   evaluation = 505,digits=0,
                   span = 1/50,digits=1,
                   family = c("symmetric", "gaussian"))

plot(LOS1,type='l',col='red',ylim=c(0,1),lwd = 2)

LOS1<-loess.smooth(seq(1,length(DivWOR)),DivWOR,degree = 1,
                   evaluation = 505,digits=0,
                   span = 1/50,digits=1,
                   family = c("symmetric", "gaussian"))
par(new=T)
plot(LOS1,type='l',col='dodgerblue',ylim=c(0,1),lwd = 2)

LOS1<-loess.smooth(seq(1,length(DivSEV)),DivSEV,degree = 1,
                   evaluation = 505,digits=0,
                   span = 1/50,digits=1,
                   family = c("symmetric", "gaussian"))
par(new=T)
plot(LOS1,type='l',col='orange',ylim=c(0,1),lwd = 2)

LOS1<-loess.smooth(seq(1,length(DivCO)),DivCO,degree = 1,
                   evaluation = 505,digits=0,
                   span = 1/50,digits=1,
                   family = c("symmetric", "gaussian"))
par(new=T)
plot(LOS1,type='l',col='green',ylim=c(0,1),lwd = 2)




AFV
LOS1<-loess.smooth(seq(1,length(AFV[,1])),AFV[,1],degree = 1,
                   evaluation = 1800,digits=0,
                   span = 1/50,digits=1,
                   family = c("symmetric", "gaussian"))

plot(LOS1,type='l',col='red',ylim=c(0,1),lwd = 2)

LOS1<-loess.smooth(seq(1,length(AFV[,1])),AFV[,2],degree = 1,
                   evaluation = 500,digits=0,
                   span = 1/100,digits=1,
                   family = c("symmetric", "gaussian"))

par(new=T)
plot(LOS1,type='l',col='dodgerblue',ylim=c(0,1),lwd = 2)

LOS1<-loess.smooth(seq(1,length(AFV[,1])),AFV[,3],degree = 1,
                   evaluation = 500,digits=0,
                   span = 1/100,digits=1,
                   family = c("symmetric", "gaussian"))

par(new=T)
plot(LOS1,type='l',col='orange',ylim=c(0,1),lwd = 2)

LOS1<-loess.smooth(seq(1,length(AFV[,1])),AFV[,4],degree = 1,
                   evaluation = 500,digits=0,
                   span = 1/100,digits=1,
                   family = c("symmetric", "gaussian"))

par(new=T)
plot(LOS1,type='l',col='green',ylim=c(0,1),lwd = 2)


#DIVgence ORFs
setwd('z:/2020/Raj_Paper/CleanData/blast_orfs/')

#DIVgence ORFs
{
  
  
  
  ORFs<-read.table('out.cds.blast')
  ORFF<-ORFs[ORFs[,1]=='giCAL',c(13,14)]
  
  lorf<-c()
  for (i in seq(1,length(ORFF[,1]))){
    print(max(c(ORFF[i,1],ORFF[i,2]))-min(c(ORFF[i,1],ORFF[i,2])))
    lorf<-c(lorf,max(c(ORFF[i,1],ORFF[i,2]))-min(c(ORFF[i,1],ORFF[i,2])))
  }
  olorf<-order(lorf)
  oolorf<-olorf[seq(1,length(lorf)-7)]
  ORFF<-ORFF[oolorf,]
  
  DivCAL<-rowMeans(Cor1[,c(1,2,3)])
  lorf<-c()
  for (p in seq(1,length(ORFF[,1]))){
    ll<-(length(seq(ORFF[p,1],ORFF[p,2])))
    lorf<-c(lorf,ll)
    
  }
  lorf
  lorfo<-order(lorf,decreasing = F)
  
  Pdf<-c(rep(0,3175))
  Pfd1<-c(0,0,0)
  for (p in lorfo){
    print(ORFF[p,])
    P1<-c(rep(0,3175))
    P1[seq(ORFF[p,1],ORFF[p,2])]=1
    Pdf<-rbind(Pdf,P1)
    print((colSums(Pdf[,seq(ORFF[p,1],ORFF[p,2])])))
    MCS<-max(colSums(Pdf[,seq(ORFF[p,1],ORFF[p,2])]))
    ORFFF<-(cbind(ORFF[p,],MCS))
    Pfd1<-rbind(Pfd1,ORFFF)
  }
  
  Pfd1<-Pfd1[-1,]
  inens<-c()
  for (p in seq(1,length(Pfd1[,1]))){
    inens<-c(inens,mean(DivCAL[seq(Pfd1[p,1],Pfd1[p,2])]))
  }
  inens<-inens*1.5
  
  #plot(NULL,ylim=c(0,40),xlim = c(0,3175))
  par(oma=c(3,0.5,3,0.5),mar=c(1,3,3,0.5))
  par(bty = 'n') 
  plot(NULL,ylim=c(-2,40),xlim = c(0,3175))
  
  
  for (p in seq(1,length(Pfd1[,1]))){
    color_transparent <- adjustcolor('red', alpha.f = inens[p])
    xbottom <- Pfd1[p,1]
    ybottom <- Pfd1[p,3]+35
    xtop <- Pfd1[p,2]
    ytop<- Pfd1[p,3]+36
    xpol <- c(xbottom,xbottom,xtop,xtop)
    ypol <- c(ybottom,ytop,ytop,ybottom)
    polygon(x=cuts(t(xpol)), y=cuts(t(ypol)), col = color_transparent, border ='black')
  }
  
  
  ORFF<-ORFs[ORFs[,1]=='giWOR',c(13,14)]
  lorf<-c()
  for (i in seq(1,length(ORFF[,1]))){
    print(max(c(ORFF[i,1],ORFF[i,2]))-min(c(ORFF[i,1],ORFF[i,2])))
    lorf<-c(lorf,max(c(ORFF[i,1],ORFF[i,2]))-min(c(ORFF[i,1],ORFF[i,2])))
  }
  olorf<-order(lorf)
  oolorf<-olorf[seq(1,length(lorf)-7)]
  ORFF<-ORFF[oolorf,]
  
  lorf<-c()
  for (p in seq(1,length(ORFF[,1]))){
    ll<-(length(seq(ORFF[p,1],ORFF[p,2])))
    lorf<-c(lorf,ll)
    
  }
  lorf
  lorfo<-order(lorf,decreasing = F)
  
  Pdf<-c(rep(0,3175))
  Pfd1<-c(0,0,0)
  for (p in lorfo){
    print(ORFF[p,])
    P1<-c(rep(0,3175))
    P1[seq(ORFF[p,1],ORFF[p,2])]=1
    Pdf<-rbind(Pdf,P1)
    print((colSums(Pdf[,seq(ORFF[p,1],ORFF[p,2])])))
    MCS<-max(colSums(Pdf[,seq(ORFF[p,1],ORFF[p,2])]))
    ORFFF<-(cbind(ORFF[p,],MCS))
    Pfd1<-rbind(Pfd1,ORFFF)
  }
  
  Pfd1<-Pfd1[-1,]
  inens<-c()
  for (p in seq(1,length(Pfd1[,1]))){
    inens<-c(inens,mean(DivCAL[seq(Pfd1[p,1],Pfd1[p,2])]))
  }
  inens<-inens*1.5
  #plot(NULL,ylim=c(0,40),xlim = c(0,3175))
  
  for (p in seq(1,length(Pfd1[,1]))){
    color_transparent <- adjustcolor('dodgerblue', alpha.f = inens[p])
    xbottom <- Pfd1[p,1]
    ybottom <- Pfd1[p,3]+27
    xtop <- Pfd1[p,2]
    ytop<- Pfd1[p,3]+28
    xpol <- c(xbottom,xbottom,xtop,xtop)
    ypol <- c(ybottom,ytop,ytop,ybottom)
    polygon(x=cuts(t(xpol)), y=cuts(t(ypol)), col = color_transparent, border ='black')
  }
  
  
  
  ORFF<-ORFs[ORFs[,1]=='giSEV',c(13,14)]
  lorf<-c()
  for (i in seq(1,length(ORFF[,1]))){
    print(max(c(ORFF[i,1],ORFF[i,2]))-min(c(ORFF[i,1],ORFF[i,2])))
    lorf<-c(lorf,max(c(ORFF[i,1],ORFF[i,2]))-min(c(ORFF[i,1],ORFF[i,2])))
  }
  olorf<-order(lorf)
  oolorf<-olorf[seq(1,length(lorf)-7)]
  ORFF<-ORFF[oolorf,]
  
  
  
  lorf<-c()
  for (p in seq(1,length(ORFs[,1]))){
    ll<-length(seq(ORFs[p,13],ORFs[p,14]))
    lorf<-c(lorf,ll)
  }
  
  lorf<-c()
  for (p in seq(1,length(ORFF[,1]))){
    ll<-(length(seq(ORFF[p,1],ORFF[p,2])))
    lorf<-c(lorf,ll)
    
  }
  lorf
  lorfo<-order(lorf,decreasing = F)
  
  Pdf<-c(rep(0,3175))
  Pfd1<-c(0,0,0)
  for (p in lorfo){
    print(ORFF[p,])
    P1<-c(rep(0,3175))
    P1[seq(ORFF[p,1],ORFF[p,2])]=1
    Pdf<-rbind(Pdf,P1)
    print((colSums(Pdf[,seq(ORFF[p,1],ORFF[p,2])])))
    MCS<-max(colSums(Pdf[,seq(ORFF[p,1],ORFF[p,2])]))
    print(MCS)
    ORFFF<-(cbind(ORFF[p,],MCS))
    Pfd1<-rbind(Pfd1,ORFFF)
  }
  
  Pfd1<-Pfd1[-1,]
  inens<-c()
  for (p in seq(1,length(Pfd1[,1]))){
    inens<-c(inens,mean(DivCAL[seq(Pfd1[p,1],Pfd1[p,2])]))
  }
  inens<-inens*1.5
  
  for (p in seq(1,length(Pfd1[,1]))){
    color_transparent <- adjustcolor('orange', alpha.f = inens[p])
    xbottom <- Pfd1[p,1]
    ybottom <- Pfd1[p,3]+19
    xtop <- Pfd1[p,2]
    ytop<- Pfd1[p,3]+20
    xpol <- c(xbottom,xbottom,xtop,xtop)
    ypol <- c(ybottom,ytop,ytop,ybottom)
    polygon(x=cuts(t(xpol)), y=cuts(t(ypol)), col = color_transparent, border ='black')
  }
  
  
  
  
  
  ORFF<-ORFs[ORFs[,1]=='giCO',c(13,14)]
  lorf<-c()
  for (i in seq(1,length(ORFF[,1]))){
    print(max(c(ORFF[i,1],ORFF[i,2]))-min(c(ORFF[i,1],ORFF[i,2])))
    lorf<-c(lorf,max(c(ORFF[i,1],ORFF[i,2]))-min(c(ORFF[i,1],ORFF[i,2])))
  }
  olorf<-order(lorf)
  oolorf<-olorf[seq(1,length(lorf)-7)]
  ORFF<-ORFF[oolorf,]
  
  lorf<-c()
  for (p in seq(1,length(ORFF[,1]))){
    ll<-(length(seq(ORFF[p,1],ORFF[p,2])))
    lorf<-c(lorf,ll)
    
  }
  lorf
  lorfo<-order(lorf,decreasing = F)
  
  Pdf<-c(rep(0,3175))
  Pfd1<-c(0,0,0)
  for (p in lorfo){
    print(ORFF[p,])
    P1<-c(rep(0,3175))
    P1[seq(ORFF[p,1],ORFF[p,2])]=1
    Pdf<-rbind(Pdf,P1)
    print((colSums(Pdf[,seq(ORFF[p,1],ORFF[p,2])])))
    MCS<-max(colSums(Pdf[,seq(ORFF[p,1],ORFF[p,2])]))
    ORFFF<-(cbind(ORFF[p,],MCS))
    Pfd1<-rbind(Pfd1,ORFFF)
  }
  
  Pfd1<-Pfd1[-1,]
  inens<-c()
  for (p in seq(1,length(Pfd1[,1]))){
    inens<-c(inens,mean(DivCAL[seq(Pfd1[p,1],Pfd1[p,2])]))
  }
  inens<-inens*1.5
  #plot(NULL,ylim=c(0,40),xlim = c(0,3175))
  
  for (p in seq(1,length(Pfd1[,1]))){
    color_transparent <- adjustcolor('green', alpha.f = inens[p])
    xbottom <- Pfd1[p,1]
    ybottom <- Pfd1[p,3]+11
    xtop <- Pfd1[p,2]
    ytop<- Pfd1[p,3]+12
    xpol <- c(xbottom,xbottom,xtop,xtop)
    ypol <- c(ybottom,ytop,ytop,ybottom)
    polygon(x=cuts(t(xpol)), y=cuts(t(ypol)), col = color_transparent, border ='black')
  }
}  

#DIVgence SNCs
{
  #plot(NULL,ylim=c(0,40),xlim = c(0,3175))
  sncRNA<-read.table('out.snc.blast')
  sncRNA<-sncRNA[sncRNA[,1]%in%sncRNACA[,2]==T,]
  Pfd1<-sncRNA[,c(9,10)]
  inens<-c()
  for (p in seq(1,length(Pfd1[,1]))){
    inens<-c(inens,mean(DivCAL[seq(Pfd1[p,1],Pfd1[p,2])]))
  }
  for (p in seq(1,length(Pfd1[,1]))){
    color_transparent <- adjustcolor('black', alpha.f = inens[p])
    xbottom <- Pfd1[p,1]
    ybottom <- 35
    xtop <- Pfd1[p,2]
    ytop<- 36
    xpol <- c(xbottom,xbottom,xtop,xtop)
    ypol <- c(ybottom,ytop,ytop,ybottom)
    polygon(x=cuts(t(xpol)), y=cuts(t(ypol)), col = color_transparent)  
  }
  sncRNA<-read.table('out.snc.blast')
  sncRNA<-sncRNA[sncRNA[,1]%in%sncRNATWOR[,2]==T,]
  Pfd1<-sncRNA[,c(9,10)]
  inens<-c()
  for (p in seq(1,length(Pfd1[,1]))){
    inens<-c(inens,mean(DivWOR[seq(Pfd1[p,1],Pfd1[p,2])]))
  }
  for (p in seq(1,length(Pfd1[,1]))){
    color_transparent <- adjustcolor('black', alpha.f = inens[p])
    xbottom <- Pfd1[p,1]
    ybottom <- 27
    xtop <- Pfd1[p,2]
    ytop<- 28
    xpol <- c(xbottom,xbottom,xtop,xtop)
    ypol <- c(ybottom,ytop,ytop,ybottom)
    polygon(x=cuts(t(xpol)), y=cuts(t(ypol)), col = color_transparent)  
  }  
  sncRNA<-read.table('out.snc.blast')
  sncRNA<-sncRNA[sncRNA[,1]%in%sncRNAS[,2]==T,]
  Pfd1<-sncRNA[,c(9,10)]
  inens<-c()
  for (p in seq(1,length(Pfd1[,1]))){
    inens<-c(inens,mean(DivSEV[seq(Pfd1[p,1],Pfd1[p,2])]))
  }
  for (p in seq(1,length(Pfd1[,1]))){
    color_transparent <- adjustcolor('black', alpha.f = inens[p])
    xbottom <- Pfd1[p,1]
    ybottom <- 19
    xtop <- Pfd1[p,2]
    ytop<- 20
    xpol <- c(xbottom,xbottom,xtop,xtop)
    ypol <- c(ybottom,ytop,ytop,ybottom)
    polygon(x=cuts(t(xpol)), y=cuts(t(ypol)), col = color_transparent)  
  }
  sncRNA<-read.table('out.snc.blast')
  sncRNA<-sncRNA[sncRNA[,1]%in%sncRNATCO[,2]==T,]
  Pfd1<-sncRNA[,c(9,10)]
  inens<-c()
  for (p in seq(1,length(Pfd1[,1]))){
    inens<-c(inens,mean(DivCO[seq(Pfd1[p,1],Pfd1[p,2])]))
  }
  for (p in seq(1,length(Pfd1[,1]))){
    color_transparent <- adjustcolor('black', alpha.f = inens[p])
    xbottom <- Pfd1[p,1]
    ybottom <- 11
    xtop <- Pfd1[p,2]
    ytop<- 12
    xpol <- c(xbottom,xbottom,xtop,xtop)
    ypol <- c(ybottom,ytop,ytop,ybottom)
    polygon(x=cuts(t(xpol)), y=cuts(t(ypol)), col = color_transparent)  
  }  
  
  
  
  
  
  #plot(NULL,ylim=c(0,40),xlim = c(0,3175))
  sncRNA<-read.table('out.snc.blast')
  sncRNA<-sncRNA[sncRNA[,1]%in%sncRNACA[,2]==T,]
  Pfd1<-sncRNA[,c(9,10)]
  inens<-c()
  for (p in seq(1,length(Pfd1[,1]))){
    inens<-c(inens,mean(DivCAL[seq(Pfd1[p,1],Pfd1[p,2])]))
  }
  for (p in seq(1,length(Pfd1[,1]))){
    color_transparent <- adjustcolor('black', alpha.f = inens[p])
    xbottom <- Pfd1[p,1]
    ybottom <- 35
    xtop <- Pfd1[p,2]
    ytop<- 36
    xpol <- c(xbottom,xbottom,xtop,xtop)
    ypol <- c(ybottom,ytop,ytop,ybottom)
    polygon(x=cuts(t(xpol)), y=cuts(t(ypol)), col = color_transparent)  
  }
  sncRNA<-read.table('out.snc.blast')
  sncRNA<-sncRNA[sncRNA[,1]%in%sncRNATWOR[,2]==T,]
  Pfd1<-sncRNA[,c(9,10)]
  inens<-c()
  for (p in seq(1,length(Pfd1[,1]))){
    inens<-c(inens,mean(DivWOR[seq(Pfd1[p,1],Pfd1[p,2])]))
  }
  for (p in seq(1,length(Pfd1[,1]))){
    color_transparent <- adjustcolor('black', alpha.f = inens[p])
    xbottom <- Pfd1[p,1]
    ybottom <- 27
    xtop <- Pfd1[p,2]
    ytop<- 28
    xpol <- c(xbottom,xbottom,xtop,xtop)
    ypol <- c(ybottom,ytop,ytop,ybottom)
    polygon(x=cuts(t(xpol)), y=cuts(t(ypol)), col = color_transparent)  
  }  
  sncRNA<-read.table('out.snc.blast')
  sncRNA<-sncRNA[sncRNA[,1]%in%sncRNAS[,2]==T,]
  Pfd1<-sncRNA[,c(9,10)]
  inens<-c()
  for (p in seq(1,length(Pfd1[,1]))){
    inens<-c(inens,mean(DivSEV[seq(Pfd1[p,1],Pfd1[p,2])]))
  }
  for (p in seq(1,length(Pfd1[,1]))){
    color_transparent <- adjustcolor('black', alpha.f = inens[p])
    xbottom <- Pfd1[p,1]
    ybottom <- 19
    xtop <- Pfd1[p,2]
    ytop<- 20
    xpol <- c(xbottom,xbottom,xtop,xtop)
    ypol <- c(ybottom,ytop,ytop,ybottom)
    polygon(x=cuts(t(xpol)), y=cuts(t(ypol)), col = color_transparent)  
  }
  sncRNA<-read.table('out.snc.blast')
  sncRNA<-sncRNA[sncRNA[,1]%in%sncRNATCO[,2]==T,]
  Pfd1<-sncRNA[,c(9,10)]
  inens<-c()
  for (p in seq(1,length(Pfd1[,1]))){
    inens<-c(inens,mean(DivCO[seq(Pfd1[p,1],Pfd1[p,2])]))
  }
  for (p in seq(1,length(Pfd1[,1]))){
    color_transparent <- adjustcolor('black', alpha.f = inens[p])
    xbottom <- Pfd1[p,1]
    ybottom <- 11
    xtop <- Pfd1[p,2]
    ytop<- 12
    xpol <- c(xbottom,xbottom,xtop,xtop)
    ypol <- c(ybottom,ytop,ytop,ybottom)
    polygon(x=cuts(t(xpol)), y=cuts(t(ypol)), col = color_transparent)  
  }  
  
  
  #plot(NULL,ylim=c(0,40),xlim = c(0,3175))
  sncRNA<-read.table('out.snc.blast')
  sncRNA<-sncRNA[sncRNA[,1]%in%sncRNACA[,2]==T,]
  Pfd1<-sncRNA[,c(9,10)]
  inens<-c()
  for (p in seq(1,length(Pfd1[,1]))){
    inens<-c(inens,mean(DivCAL[seq(Pfd1[p,1],Pfd1[p,2])]))
  }
  for (p in seq(1,length(Pfd1[,1]))){
    color_transparent <- adjustcolor('black', alpha.f = .1)
    xbottom <- Pfd1[p,1]
    ybottom <- 35
    xtop <- Pfd1[p,2]
    ytop<- 36
    xpol <- c(xbottom,xbottom,xtop,xtop)
    ypol <- c(ybottom,ytop,ytop,ybottom)
    polygon(x=cuts(t(xpol)), y=cuts(t(ypol)), col = color_transparent)  
  }
  sncRNA<-read.table('out.snc.blast')
  sncRNA<-sncRNA[sncRNA[,1]%in%sncRNATWOR[,2]==T,]
  Pfd1<-sncRNA[,c(9,10)]
  inens<-c()
  for (p in seq(1,length(Pfd1[,1]))){
    inens<-c(inens,mean(DivWOR[seq(Pfd1[p,1],Pfd1[p,2])]))
  }
  for (p in seq(1,length(Pfd1[,1]))){
    color_transparent <- adjustcolor('black', alpha.f = .1)
    xbottom <- Pfd1[p,1]
    ybottom <- 27
    xtop <- Pfd1[p,2]
    ytop<- 28
    xpol <- c(xbottom,xbottom,xtop,xtop)
    ypol <- c(ybottom,ytop,ytop,ybottom)
    polygon(x=cuts(t(xpol)), y=cuts(t(ypol)), col = color_transparent)  
  }  
  sncRNA<-read.table('out.snc.blast')
  sncRNA<-sncRNA[sncRNA[,1]%in%sncRNAS[,2]==T,]
  Pfd1<-sncRNA[,c(9,10)]
  inens<-c()
  for (p in seq(1,length(Pfd1[,1]))){
    inens<-c(inens,mean(DivSEV[seq(Pfd1[p,1],Pfd1[p,2])]))
  }
  for (p in seq(1,length(Pfd1[,1]))){
    color_transparent <- adjustcolor('black', alpha.f = .1)
    xbottom <- Pfd1[p,1]
    ybottom <- 19
    xtop <- Pfd1[p,2]
    ytop<- 20
    xpol <- c(xbottom,xbottom,xtop,xtop)
    ypol <- c(ybottom,ytop,ytop,ybottom)
    polygon(x=cuts(t(xpol)), y=cuts(t(ypol)), col = color_transparent)  
  }
  sncRNA<-read.table('out.snc.blast')
  sncRNA<-sncRNA[sncRNA[,1]%in%sncRNATCO[,2]==T,]
  Pfd1<-sncRNA[,c(9,10)]
  inens<-c()
  for (p in seq(1,length(Pfd1[,1]))){
    inens<-c(inens,mean(DivCO[seq(Pfd1[p,1],Pfd1[p,2])]))
  }
  for (p in seq(1,length(Pfd1[,1]))){
    color_transparent <- adjustcolor('black', alpha.f = .1)
    xbottom <- Pfd1[p,1]
    ybottom <- 11
    xtop <- Pfd1[p,2]
    ytop<- 12
    xpol <- c(xbottom,xbottom,xtop,xtop)
    ypol <- c(ybottom,ytop,ytop,ybottom)
    polygon(x=cuts(t(xpol)), y=cuts(t(ypol)), col = color_transparent)  
  }  
  
}


#BCTV genes
{
  color_transparent <- adjustcolor('grey0', alpha.f = 0.5)
  sncRNA<-read.table('BCTVgenes_prot.blast.out')
  Gm<-c(min(sncRNA[sncRNA[,1]=='V1',c(9,10)]))
  GM<-c(max(sncRNA[sncRNA[,1]=='V1',c(9,10)]))
  xbottom <- Gm
  ybottom <- 5
  xtop <- GM
  ytop<- 6
  xpol <- c(xbottom,xbottom,xtop,xtop)
  ypol <- c(ybottom,ytop,ytop,ybottom)
  polygon(x=cuts(t(xpol)), y=cuts(t(ypol)), col = color_transparent, border ='black')
  
  Gm<-c(min(sncRNA[sncRNA[,1]=='V2',c(9,10)]))
  GM<-c(max(sncRNA[sncRNA[,1]=='V2',c(9,10)]))
  xbottom <- Gm
  ybottom <- 6
  xtop <- GM
  ytop<- 7
  xpol <- c(xbottom,xbottom,xtop,xtop)
  ypol <- c(ybottom,ytop,ytop,ybottom)
  polygon(x=cuts(t(xpol)), y=cuts(t(ypol)), col = color_transparent, border ='black')
  
  Gm<-c(min(sncRNA[sncRNA[,1]=='V3',c(9,10)]))
  GM<-c(max(sncRNA[sncRNA[,1]=='V3',c(9,10)]))
  xbottom <- Gm
  ybottom <- 5
  xtop <- GM
  ytop<- 6
  xpol <- c(xbottom,xbottom,xtop,xtop)
  ypol <- c(ybottom,ytop,ytop,ybottom)
  polygon(x=cuts(t(xpol)), y=cuts(t(ypol)), col = color_transparent, border ='black')
  
  color_transparent <- adjustcolor('darkslategray', alpha.f = 0.5)
  Gm<-c(min(sncRNA[sncRNA[,1]=='C1',c(9,10)]))
  GM<-c(max(sncRNA[sncRNA[,1]=='C1',c(9,10)]))
  xbottom <- Gm
  ybottom <- 5
  xtop <- GM
  ytop<- 6
  xpol <- c(xbottom,xbottom,xtop,xtop)
  ypol <- c(ybottom,ytop,ytop,ybottom)
  polygon(x=cuts(t(xpol)), y=cuts(t(ypol)), col = color_transparent, border ='black')
  
  Gm<-c(min(sncRNA[sncRNA[,1]=='C2',c(9,10)]))
  GM<-c(max(sncRNA[sncRNA[,1]=='C2',c(9,10)]))
  xbottom <- Gm
  ybottom <- 6
  xtop <- GM
  ytop<- 7
  xpol <- c(xbottom,xbottom,xtop,xtop)
  ypol <- c(ybottom,ytop,ytop,ybottom)
  polygon(x=cuts(t(xpol)), y=cuts(t(ypol)), col = color_transparent, border ='black')
  
  sncRNA[sncRNA[,1]=='C3',c(9,10)]
  Gm<-c(min(sncRNA[sncRNA[,1]=='C3',c(9,10)]))
  GM<-c(max(sncRNA[sncRNA[,1]=='C3',c(9,10)]))
  xbottom <- Gm
  ybottom <- 5
  xtop <- GM
  ytop<- 6
  xpol <- c(xbottom,xbottom,xtop,xtop)
  ypol <- c(ybottom,ytop,ytop,ybottom)
  polygon(x=cuts(t(xpol)), y=cuts(t(ypol)), col = color_transparent, border ='black')
  
  sncRNA[sncRNA[,1]=='C4',c(9,10)]
  Gm<-c(min(sncRNA[sncRNA[,1]=='C4',c(9,10)]))
  GM<-c(max(sncRNA[sncRNA[,1]=='C4',c(9,10)]))
  xbottom <- Gm
  ybottom <- 6
  xtop <- GM
  ytop<- 7
  xpol <- c(xbottom,xbottom,xtop,xtop)
  ypol <- c(ybottom,ytop,ytop,ybottom)
  polygon(x=cuts(t(xpol)), y=cuts(t(ypol)), col = color_transparent, border ='black')
}


#OFDdivtables
{
  
  ORFs<-read.table('out.cds.blast')
  ORFF<-ORFs[ORFs[,1]=='giCAL',c(5,13,14)]
  DivCAL<-rowMeans(Cor1[,c(1,2,3)])
  pdf<-c(0,0,0,0,0,0,0)
  for (i in seq(1,length(ORFF[,1]))){
    pp<-mean(DivCAL[seq(ORFF[i,2],ORFF[i,3])])
    ssd<-sd(DivCAL[seq(ORFF[i,2],ORFF[i,3])])
    m<-min(DivCAL[seq(ORFF[i,2],ORFF[i,3])])
    M<-max(DivCAL[seq(ORFF[i,2],ORFF[i,3])])
    pdf<-rbind(pdf,cbind('CAL',ORFF[i,],pp,ssd,m,M))
  }
  
  pdf<-pdf[-1,]
  pdf
  pdf1<-pdf
  
  ORFF<-ORFs[ORFs[,1]=='giWOR',c(5,13,14)]
  DivWOR<-rowMeans(Cor1[,c(4,5,6)])
  pdf<-c(0,0,0,0,0,0,0)
  for (i in seq(1,length(ORFF[,1]))){
    pp<-mean(DivWOR[seq(ORFF[i,2],ORFF[i,3])])
    ssd<-sd(DivWOR[seq(ORFF[i,2],ORFF[i,3])])
    m<-min(DivWOR[seq(ORFF[i,2],ORFF[i,3])])
    M<-max(DivWOR[seq(ORFF[i,2],ORFF[i,3])])
    pdf<-rbind(pdf,cbind('WOR',ORFF[i,],pp,ssd,m,M))
  }
  
  pdf<-pdf[-1,]
  pdf
  
  
  
  ORFF<-ORFs[ORFs[,1]=='giSEV',c(5,13,14)]
  DivSEV<-rowMeans(Cor1[,c(7,8,9)])
  pdf<-c(0,0,0,0,0,0,0)
  for (i in seq(1,length(ORFF[,1]))){
    pp<-mean(DivSEV[seq(ORFF[i,2],ORFF[i,3])])
    ssd<-sd(DivSEV[seq(ORFF[i,2],ORFF[i,3])])
    m<-min(DivSEV[seq(ORFF[i,2],ORFF[i,3])])
    M<-max(DivSEV[seq(ORFF[i,2],ORFF[i,3])])
    pdf<-rbind(pdf,cbind('SEV',ORFF[i,],pp,ssd,m,M))
  }
  
  pdf<-pdf[-1,]
  pdf
  
  
  ORFF<-ORFs[ORFs[,1]=='giCO',c(5,13,14)]
  DivCO<-rowMeans(Cor1[,c(10,11,12)])
  pdf<-c(0,0,0,0,0,0,0)
  for (i in seq(1,length(ORFF[,1]))){
    pp<-mean(DivCO[seq(ORFF[i,2],ORFF[i,3])])
    ssd<-sd(DivCO[seq(ORFF[i,2],ORFF[i,3])])
    m<-min(DivCO[seq(ORFF[i,2],ORFF[i,3])])
    M<-max(DivCO[seq(ORFF[i,2],ORFF[i,3])])
    pdf<-rbind(pdf,cbind('CO',ORFF[i,],pp,ssd,m,M))
  }
  
  pdf<-pdf[-1,]
  pdf
  
}

#SNCdivtables
{
  sncRNA<-read.table('out.snc.blast')
  sncRNA<-sncRNA[sncRNA[,1]%in%sncRNACA[,2]==T,]
  Pfd1<-sncRNA[,c(1,9,10)]
  inens<-c()
  inens1<-c()
  inens2<-c()
  inens3<-c()
  for (p in seq(1,length(Pfd1[,1]))){
    inens<-c(inens,mean(DivCAL[seq(Pfd1[p,2],Pfd1[p,3])]))
    inens1<-c(inens1,sd(DivCAL[seq(Pfd1[p,2],Pfd1[p,3])]))
    inens2<-c(inens2,min(DivCAL[seq(Pfd1[p,2],Pfd1[p,3])]))
    inens3<-c(inens3,max(DivCAL[seq(Pfd1[p,2],Pfd1[p,3])]))
  }
  
  cbind(Pfd1,inens,inens1,inens2,inens3)
  
  sncRNA<-read.table('out.snc.blast')
  sncRNA<-sncRNA[sncRNA[,1]%in%sncRNATWOR[,2]==T,]
  Pfd1<-sncRNA[,c(1,9,10)]
  inens<-c()
  inens1<-c()
  inens2<-c()
  inens3<-c()
  for (p in seq(1,length(Pfd1[,1]))){
    inens<-c(inens,mean(DivWOR[seq(Pfd1[p,2],Pfd1[p,3])]))
    inens1<-c(inens1,sd(DivWOR[seq(Pfd1[p,2],Pfd1[p,3])]))
    inens2<-c(inens2,min(DivWOR[seq(Pfd1[p,2],Pfd1[p,3])]))
    inens3<-c(inens3,max(DivWOR[seq(Pfd1[p,2],Pfd1[p,3])]))
  }
  
  cbind(Pfd1,inens,inens1,inens2,inens3)
  
  
  sncRNA<-read.table('out.snc.blast')
  sncRNA<-sncRNA[sncRNA[,1]%in%sncRNAS[,2]==T,]
  Pfd1<-sncRNA[,c(1,9,10)]
  inens<-c()
  inens1<-c()
  inens2<-c()
  inens3<-c()
  for (p in seq(1,length(Pfd1[,1]))){
    inens<-c(inens,mean(DivSEV[seq(Pfd1[p,2],Pfd1[p,3])]))
    inens1<-c(inens1,sd(DivSEV[seq(Pfd1[p,2],Pfd1[p,3])]))
    inens2<-c(inens2,min(DivSEV[seq(Pfd1[p,2],Pfd1[p,3])]))
    inens3<-c(inens3,max(DivSEV[seq(Pfd1[p,2],Pfd1[p,3])]))
  }
  
  cbind(Pfd1,inens,inens1,inens2,inens3)
  
  
  sncRNA<-read.table('out.snc.blast')
  sncRNA<-sncRNA[sncRNA[,1]%in%sncRNATCO[,2]==T,]
  Pfd1<-sncRNA[,c(1,9,10)]
  inens<-c()
  inens1<-c()
  inens2<-c()
  inens3<-c()
  for (p in seq(1,length(Pfd1[,1]))){
    inens<-c(inens,mean(DivCO[seq(Pfd1[p,2],Pfd1[p,3])]))
    inens1<-c(inens1,sd(DivCO[seq(Pfd1[p,2],Pfd1[p,3])]))
    inens2<-c(inens2,min(DivCO[seq(Pfd1[p,2],Pfd1[p,3])]))
    inens3<-c(inens3,max(DivCO[seq(Pfd1[p,2],Pfd1[p,3])]))
  }
  
  cbind(Pfd1,inens,inens1,inens2,inens3)
  
}

#Genedivtables
{
  sncRNA<-read.table('BCTVgenes_prot.blast.out')
  genes<-c(0,0,0)
  for (i in unique(sncRNA[,1])){
    g1<-(cbind(i,min(sncRNA[sncRNA[,1]==i,c(9,10)]),max(sncRNA[sncRNA[,1]==i,c(9,10)])))
    genes<-rbind(genes,g1)
  }
  genes
  genes<-genes[-1,]
  inens<-c()
  inens1<-c()
  inens2<-c()
  inens3<-c()
  
  for (p in seq(1,length(genes[,1]))){
    inens<-c(inens,mean(DivCAL[seq(genes[p,2],genes[p,3])]))
    inens1<-c(inens1,sd(DivCAL[seq(genes[p,2],genes[p,3])]))
    inens2<-c(inens2,min(DivCAL[seq(genes[p,2],genes[p,3])]))
    inens3<-c(inens3,max(DivCAL[seq(genes[p,2],genes[p,3])]))
  }
  as.data.frame(cbind(genes,inens,inens1,inens2,inens3))
  
  
  sncRNA<-read.table('BCTVgenes_prot.blast.out')
  genes<-c(0,0,0)
  for (i in unique(sncRNA[,1])){
    g1<-(cbind(i,min(sncRNA[sncRNA[,1]==i,c(9,10)]),max(sncRNA[sncRNA[,1]==i,c(9,10)])))
    genes<-rbind(genes,g1)
  }
  genes
  genes<-genes[-1,]
  inens<-c()
  inens1<-c()
  inens2<-c()
  inens3<-c()
  
  for (p in seq(1,length(genes[,1]))){
    inens<-c(inens,mean(DivWOR[seq(genes[p,2],genes[p,3])]))
    inens1<-c(inens1,sd(DivWOR[seq(genes[p,2],genes[p,3])]))
    inens2<-c(inens2,min(DivWOR[seq(genes[p,2],genes[p,3])]))
    inens3<-c(inens3,max(DivWOR[seq(genes[p,2],genes[p,3])]))
  }
  as.data.frame(cbind(genes,inens,inens1,inens2,inens3))  
  
  sncRNA<-read.table('BCTVgenes_prot.blast.out')
  genes<-c(0,0,0)
  for (i in unique(sncRNA[,1])){
    g1<-(cbind(i,min(sncRNA[sncRNA[,1]==i,c(9,10)]),max(sncRNA[sncRNA[,1]==i,c(9,10)])))
    genes<-rbind(genes,g1)
  }
  genes
  genes<-genes[-1,]
  inens<-c()
  inens1<-c()
  inens2<-c()
  inens3<-c()
  
  for (p in seq(1,length(genes[,1]))){
    inens<-c(inens,mean(DivSEV[seq(genes[p,2],genes[p,3])]))
    inens1<-c(inens1,sd(DivSEV[seq(genes[p,2],genes[p,3])]))
    inens2<-c(inens2,min(DivSEV[seq(genes[p,2],genes[p,3])]))
    inens3<-c(inens3,max(DivSEV[seq(genes[p,2],genes[p,3])]))
  }
  as.data.frame(cbind(genes,inens,inens1,inens2,inens3))
  
  
  
  sncRNA<-read.table('BCTVgenes_prot.blast.out')
  genes<-c(0,0,0)
  for (i in unique(sncRNA[,1])){
    g1<-(cbind(i,min(sncRNA[sncRNA[,1]==i,c(9,10)]),max(sncRNA[sncRNA[,1]==i,c(9,10)])))
    genes<-rbind(genes,g1)
  }
  genes
  genes<-genes[-1,]
  inens<-c()
  inens1<-c()
  inens2<-c()
  inens3<-c()
  
  for (p in seq(1,length(genes[,1]))){
    inens<-c(inens,mean(DivCO[seq(genes[p,2],genes[p,3])]))
    inens1<-c(inens1,sd(DivCO[seq(genes[p,2],genes[p,3])]))
    inens2<-c(inens2,min(DivCO[seq(genes[p,2],genes[p,3])]))
    inens3<-c(inens3,max(DivCO[seq(genes[p,2],genes[p,3])]))
  }
  as.data.frame(cbind(genes,inens,inens1,inens2,inens3))
  
}  


###HS###
DivCAL<-AFV[,1]
DivWOR<-AFV[,2]
DivSEV<-AFV[,3]
DivCO<-AFV[,4]

#OFDdivtables
{
  
  ORFs<-read.table('out.cds.blast')
  ORFF<-ORFs[ORFs[,1]=='giCAL',c(5,13,14)]
  pdf<-c(0,0,0,0,0,0,0)
  for (i in seq(1,length(ORFF[,1]))){
    pp<-mean(DivCAL[seq(ORFF[i,2],ORFF[i,3])])
    ssd<-sd(DivCAL[seq(ORFF[i,2],ORFF[i,3])])
    m<-min(DivCAL[seq(ORFF[i,2],ORFF[i,3])])
    M<-max(DivCAL[seq(ORFF[i,2],ORFF[i,3])])
    pdf<-rbind(pdf,cbind('CAL',ORFF[i,],pp,ssd,m,M))
  }
  
  pdf<-pdf[-1,]
  pdf
  pdf1<-pdf
  
  ORFF<-ORFs[ORFs[,1]=='giWOR',c(5,13,14)]
  pdf<-c(0,0,0,0,0,0,0)
  for (i in seq(1,length(ORFF[,1]))){
    pp<-mean(DivWOR[seq(ORFF[i,2],ORFF[i,3])])
    ssd<-sd(DivWOR[seq(ORFF[i,2],ORFF[i,3])])
    m<-min(DivWOR[seq(ORFF[i,2],ORFF[i,3])])
    M<-max(DivWOR[seq(ORFF[i,2],ORFF[i,3])])
    pdf<-rbind(pdf,cbind('WOR',ORFF[i,],pp,ssd,m,M))
  }
  
  pdf<-pdf[-1,]
  pdf
  
  
  
  ORFF<-ORFs[ORFs[,1]=='giSEV',c(5,13,14)]
  pdf<-c(0,0,0,0,0,0,0)
  for (i in seq(1,length(ORFF[,1]))){
    pp<-mean(DivSEV[seq(ORFF[i,2],ORFF[i,3])])
    ssd<-sd(DivSEV[seq(ORFF[i,2],ORFF[i,3])])
    m<-min(DivSEV[seq(ORFF[i,2],ORFF[i,3])])
    M<-max(DivSEV[seq(ORFF[i,2],ORFF[i,3])])
    pdf<-rbind(pdf,cbind('SEV',ORFF[i,],pp,ssd,m,M))
  }
  
  pdf<-pdf[-1,]
  pdf
  
  
  ORFF<-ORFs[ORFs[,1]=='giCO',c(5,13,14)]
  pdf<-c(0,0,0,0,0,0,0)
  for (i in seq(1,length(ORFF[,1]))){
    pp<-mean(DivCO[seq(ORFF[i,2],ORFF[i,3])])
    ssd<-sd(DivCO[seq(ORFF[i,2],ORFF[i,3])])
    m<-min(DivCO[seq(ORFF[i,2],ORFF[i,3])])
    M<-max(DivCO[seq(ORFF[i,2],ORFF[i,3])])
    pdf<-rbind(pdf,cbind('CO',ORFF[i,],pp,ssd,m,M))
  }
  
  pdf<-pdf[-1,]
  pdf
  
}

#SNCdivtables
{
  sncRNA<-read.table('out.snc.blast')
  sncRNA<-sncRNA[sncRNA[,1]%in%sncRNACA[,2]==T,]
  Pfd1<-sncRNA[,c(1,9,10)]
  inens<-c()
  inens1<-c()
  inens2<-c()
  inens3<-c()
  for (p in seq(1,length(Pfd1[,1]))){
    inens<-c(inens,mean(DivCAL[seq(Pfd1[p,2],Pfd1[p,3])]))
    inens1<-c(inens1,sd(DivCAL[seq(Pfd1[p,2],Pfd1[p,3])]))
    inens2<-c(inens2,min(DivCAL[seq(Pfd1[p,2],Pfd1[p,3])]))
    inens3<-c(inens3,max(DivCAL[seq(Pfd1[p,2],Pfd1[p,3])]))
  }
  
  cbind(Pfd1,inens,inens1,inens2,inens3)
  
  sncRNA<-read.table('out.snc.blast')
  sncRNA<-sncRNA[sncRNA[,1]%in%sncRNATWOR[,2]==T,]
  Pfd1<-sncRNA[,c(1,9,10)]
  inens<-c()
  inens1<-c()
  inens2<-c()
  inens3<-c()
  for (p in seq(1,length(Pfd1[,1]))){
    inens<-c(inens,mean(DivWOR[seq(Pfd1[p,2],Pfd1[p,3])]))
    inens1<-c(inens1,sd(DivWOR[seq(Pfd1[p,2],Pfd1[p,3])]))
    inens2<-c(inens2,min(DivWOR[seq(Pfd1[p,2],Pfd1[p,3])]))
    inens3<-c(inens3,max(DivWOR[seq(Pfd1[p,2],Pfd1[p,3])]))
  }
  
  cbind(Pfd1,inens,inens1,inens2,inens3)
  
  
  sncRNA<-read.table('out.snc.blast')
  sncRNA<-sncRNA[sncRNA[,1]%in%sncRNAS[,2]==T,]
  Pfd1<-sncRNA[,c(1,9,10)]
  inens<-c()
  inens1<-c()
  inens2<-c()
  inens3<-c()
  for (p in seq(1,length(Pfd1[,1]))){
    inens<-c(inens,mean(DivSEV[seq(Pfd1[p,2],Pfd1[p,3])]))
    inens1<-c(inens1,sd(DivSEV[seq(Pfd1[p,2],Pfd1[p,3])]))
    inens2<-c(inens2,min(DivSEV[seq(Pfd1[p,2],Pfd1[p,3])]))
    inens3<-c(inens3,max(DivSEV[seq(Pfd1[p,2],Pfd1[p,3])]))
  }
  
  cbind(Pfd1,inens,inens1,inens2,inens3)
  
  
  sncRNA<-read.table('out.snc.blast')
  sncRNA<-sncRNA[sncRNA[,1]%in%sncRNATCO[,2]==T,]
  Pfd1<-sncRNA[,c(1,9,10)]
  inens<-c()
  inens1<-c()
  inens2<-c()
  inens3<-c()
  for (p in seq(1,length(Pfd1[,1]))){
    inens<-c(inens,mean(DivCO[seq(Pfd1[p,2],Pfd1[p,3])]))
    inens1<-c(inens1,sd(DivCO[seq(Pfd1[p,2],Pfd1[p,3])]))
    inens2<-c(inens2,min(DivCO[seq(Pfd1[p,2],Pfd1[p,3])]))
    inens3<-c(inens3,max(DivCO[seq(Pfd1[p,2],Pfd1[p,3])]))
  }
  
  cbind(Pfd1,inens,inens1,inens2,inens3)
  
}

#Genedivtables
{
  sncRNA<-read.table('BCTVgenes_prot.blast.out')
  genes<-c(0,0,0)
  for (i in unique(sncRNA[,1])){
    g1<-(cbind(i,min(sncRNA[sncRNA[,1]==i,c(9,10)]),max(sncRNA[sncRNA[,1]==i,c(9,10)])))
    genes<-rbind(genes,g1)
  }
  genes
  genes<-genes[-1,]
  inens<-c()
  inens1<-c()
  inens2<-c()
  inens3<-c()
  
  for (p in seq(1,length(genes[,1]))){
    inens<-c(inens,mean(DivCAL[seq(genes[p,2],genes[p,3])]))
    inens1<-c(inens1,sd(DivCAL[seq(genes[p,2],genes[p,3])]))
    inens2<-c(inens2,min(DivCAL[seq(genes[p,2],genes[p,3])]))
    inens3<-c(inens3,max(DivCAL[seq(genes[p,2],genes[p,3])]))
  }
  as.data.frame(cbind(genes,inens,inens1,inens2,inens3))
  
  
  sncRNA<-read.table('BCTVgenes_prot.blast.out')
  genes<-c(0,0,0)
  for (i in unique(sncRNA[,1])){
    g1<-(cbind(i,min(sncRNA[sncRNA[,1]==i,c(9,10)]),max(sncRNA[sncRNA[,1]==i,c(9,10)])))
    genes<-rbind(genes,g1)
  }
  genes
  genes<-genes[-1,]
  inens<-c()
  inens1<-c()
  inens2<-c()
  inens3<-c()
  
  for (p in seq(1,length(genes[,1]))){
    inens<-c(inens,mean(DivWOR[seq(genes[p,2],genes[p,3])]))
    inens1<-c(inens1,sd(DivWOR[seq(genes[p,2],genes[p,3])]))
    inens2<-c(inens2,min(DivWOR[seq(genes[p,2],genes[p,3])]))
    inens3<-c(inens3,max(DivWOR[seq(genes[p,2],genes[p,3])]))
  }
  as.data.frame(cbind(genes,inens,inens1,inens2,inens3))  
  
  sncRNA<-read.table('BCTVgenes_prot.blast.out')
  genes<-c(0,0,0)
  for (i in unique(sncRNA[,1])){
    g1<-(cbind(i,min(sncRNA[sncRNA[,1]==i,c(9,10)]),max(sncRNA[sncRNA[,1]==i,c(9,10)])))
    genes<-rbind(genes,g1)
  }
  genes
  genes<-genes[-1,]
  inens<-c()
  inens1<-c()
  inens2<-c()
  inens3<-c()
  
  for (p in seq(1,length(genes[,1]))){
    inens<-c(inens,mean(DivSEV[seq(genes[p,2],genes[p,3])]))
    inens1<-c(inens1,sd(DivSEV[seq(genes[p,2],genes[p,3])]))
    inens2<-c(inens2,min(DivSEV[seq(genes[p,2],genes[p,3])]))
    inens3<-c(inens3,max(DivSEV[seq(genes[p,2],genes[p,3])]))
  }
  as.data.frame(cbind(genes,inens,inens1,inens2,inens3))
  
  
  
  sncRNA<-read.table('BCTVgenes_prot.blast.out')
  genes<-c(0,0,0)
  for (i in unique(sncRNA[,1])){
    g1<-(cbind(i,min(sncRNA[sncRNA[,1]==i,c(9,10)]),max(sncRNA[sncRNA[,1]==i,c(9,10)])))
    genes<-rbind(genes,g1)
  }
  genes
  genes<-genes[-1,]
  inens<-c()
  inens1<-c()
  inens2<-c()
  inens3<-c()
  
  for (p in seq(1,length(genes[,1]))){
    inens<-c(inens,mean(DivCO[seq(genes[p,2],genes[p,3])]))
    inens1<-c(inens1,sd(DivCO[seq(genes[p,2],genes[p,3])]))
    inens2<-c(inens2,min(DivCO[seq(genes[p,2],genes[p,3])]))
    inens3<-c(inens3,max(DivCO[seq(genes[p,2],genes[p,3])]))
  }
  as.data.frame(cbind(genes,inens,inens1,inens2,inens3))
  
}  













cbind(SAMP,AFV[,2])


library(ade4)
unique(DepthMD[,2])
ROWNamova<-rownames(DepthMD[DepthMD[,2]=='BCTV-CO',])
Structure<-as.data.frame(as.vector(DepthMD[DepthMD[,2]=='BCTV-CO',4]))
colnames(Structure)<-'State'
SAMP<-t(as.data.frame(Depth2[ROWNamova,]))
colnames(SAMP)<-seq(1,length(SAMP[1,]))
SAMP1<-cbind(SAMP,AFV[,4])
SAMP2<-SAMP1[SAMP1[,20]>0,]
SAMP3<-SAMP2[,-36]
SAMP4<-as.matrix(SAMP3)
SAMP5<-c()
for (i in seq(1,length(SAMP4[,1]))){
  if (length(unique(SAMP4[i,]))>1){
    SAMP5<-rbind(SAMP5,SAMP4[i,])
  }
}

SMAP5<-data.frame(SAMP5)
DIST<-dist(t(SAMP5),method='euclidean')
heatmap(SAMP5,Rowv = NA)
library("pheatmap")
pheatmap(SAMP5,cluster_rows=F,clustering_method = 'complete',breaks = 10)


for (i in seq(1,length(SAMP5[,1]))){
  AMOV<-amova(samples = SAMP5[i,],structures = Structure)
  print(AMOV$componentsofcovariance)
}

AMOV<-amova(samples=data.frame(SAMP5),structures = Structure)
AMOV$results
AMOV$componentsofcovariance
AMOV$statphi



data(humDNAm)
dpcoahum <- dpcoa(data.frame(t(humDNAm$samples)), 
                  sqrt(humDNAm$distances),scan = FALSE, nf = 2)
?dpcoa
plot(dpcoahum)





unique(as.vector(SAMP3[2,]))

DIST<-dist(t(SAMP5),method='euclidean')

dpcoahum <- dpcoa(data.frame(SAMP5),
                  sqrt(DIST), scan = FALSE, nf = 300)
plot(dpcoahum,col =as.numeric(Structure))


par(mfrow=c(1,1))
d <- dist(t(SAMP3)) # euclidean distances between the rows
fit <- cmdscale(d,eig=TRUE, k=2) # k is the number of dim
fit # view results

# plot solution
x <- fit$points[,1]
y <- fit$points[,2]

SI<-as.data.frame(unique(as.numeric(Structure[,1])))
SI
orderN<-order(SI)
color_colors<-rainbow(6)

plot(x, y, xlab="Coordinate 1", ylab="Coordinate 2",
     main="BCTV Virus Genome Sequences", type="p",pch=19,cex=1.5,col=color_colors[as.numeric(Structure[,1])])
legend(x = "topright",          # Position
       legend = unique(Structure[,1])[orderN],  # Legend texts
       pch = 19,           # Line types
       col = color_colors,           # Line colors
       lwd = 0)                 # Line width


install.packages("dendextend")
library(dendextend)
hc<-hclust(DIST)
dhc <- as.dendrogram(hc)
plot(dhc)
dend2 <- cut(dhc, h = 17)
dend2 <- color_branches(dhc, k = 12)
plot(dend2)

# creates a single item vector of the clusters    
myclusters <- cutree(dhc, k=12)

# make the dataframe of two columns cluster number and label
clusterDF <-  data.frame(Cluster = as.numeric(unlist(myclusters)),
                         Branch = names(myclusters))
clusterDF[order(clusterDF[,1]),]

# sort by cluster ascending
clusterDFSort <- clusterDF %>% arrange(Cluster)
Structure


##Hierachcal levels
##Strain
library(ade4)
unique(DepthMD[,2])
ROWNamova<-rownames(DepthMD)
Structure<-as.data.frame(as.vector(DepthMD[,2]))
colnames(Structure)<-'Strain'
SAMP<-t(as.data.frame(Depth2[ROWNamova,]))
colnames(SAMP)<-seq(1,length(SAMP[1,]))
SAMP1<-cbind(SAMP,AFV[,4])
SAMP2<-SAMP1[SAMP1[,20]>0,]
SAMP3<-SAMP2[,-112]
SAMP4<-as.matrix(SAMP3)
SAMP5<-c()
for (i in seq(1,length(SAMP4[,1]))){
  if (length(unique(SAMP4[i,]))>1){
    SAMP5<-rbind(SAMP5,SAMP4[i,])
  }
}

mean(rowMeans(Cor1[,c(1,2,3)]))
mean(rowMeans(Cor1[,c(4,5,6)]))
mean(rowMeans(Cor1[,c(7,8,9)]))
mean(rowMeans(Cor1[,c(10,11,12)]))
hist(rowMeans(Cor1[,c(1,2,3)]))
hist(rowMeans(Cor1[,c(4,5,6)]))
hist(rowMeans(Cor1[,c(7,8,9)]))
hist(rowMeans(Cor1[,c(10,11,12)]))

VectorRM<-rowMeans(Cor1[,c(1,2,3)])
VectorRM<-data.frame(VectorRM)
VectorRM<-cbind(VectorRM,VectorRM)
VectorRM[VectorRM[,1]>0.5,]
length(VectorRM[VectorRM[,1]>0.5,1])


AMOV<-amova(samples=data.frame(SAMP5),structures = Structure)
AMOV$results
AMOV$componentsofcovariance
AMOV$statphi

for (i in seq(1,length(SAMP5[,1]))){
  AMOV<-amova(samples = data.frame(t(SAMP5[i,])),structures = Structure)
  print(AMOV$componentsofcovariance[1,])
}




#GEOGRAPHIC LOCATION

#COLORADO
{
  unique(CORN[,4])
  
  FstALL<-data.frame()
  Cor1<-data.frame()
  AFV<-data.frame(1,2,3,4,5,6)
  colnames(AFV)<-c('CO','ID','NE','OR','CA','UN')
  ID<-CORN[CORN[,4]=='ID',]
  NE<-CORN[CORN[,4]=='NE',]
  OR<-CORN[CORN[,4]=='OR',]
  CA<-CORN[CORN[,4]=='CA',]
  CO<-CORN[CORN[,4]=='CO',]
  UN<-CORN[CORN[,4]=='UN',]
  
  for (i in seq(1,3175)){
    i<-i
    print(i)
    X<-as.vector(Depth2[rownames(ID),i])
    Xa<-(length(X[X==1])/length(X))^2
    Xt<-(length(X[X==2])/length(X))^2
    Xg<-(length(X[X==3])/length(X))^2
    Xc<-(length(X[X==4])/length(X))^2
    Xx<-(length(X[X==5])/length(X))^2
    SX<-1-sum(Xa,Xt,Xg,Xc,Xx)
    V1<-c(Xa,Xt,Xg,Xc,Xx)
    SUMX<-(SX)
    
    X1<-as.vector(Depth2[rownames(NE),i])
    Xa1<-(length(X1[X1==1])/length(X1))^2
    Xt1<-(length(X1[X1==2])/length(X1))^2
    Xg1<-(length(X1[X1==3])/length(X1))^2
    Xc1<-(length(X1[X1==4])/length(X1))^2
    Xx1<-(length(X1[X1==5])/length(X1))^2
    SX1<-1-sum(Xa1,Xt1,Xg1,Xc1,Xx1)
    SUMX1<-(SX1)
    FstSite<-((SX+SX1)-(SX1))/(SX+SX1)
    V2<-c(Xa1,Xt1,Xg1,Xc1,Xx1)
    C1<-cor(V1,V2)
    
    X2<-as.vector(Depth2[rownames(OR),i])
    Xa2<-(length(X2[X2==1])/length(X2))^2
    Xt2<-(length(X2[X2==2])/length(X2))^2
    Xg2<-(length(X2[X2==3])/length(X2))^2
    Xc2<-(length(X2[X2==4])/length(X2))^2
    Xx2<-(length(X2[X2==5])/length(X2))^2
    SX2<-1-sum(Xa2,Xt2,Xg2,Xc2,Xx2)
    V3<-c(Xa2,Xt2,Xg2,Xc2,Xx2)
    
    X3<-as.vector(Depth2[rownames(CA),i])
    Xa3<-(length(X3[X3==1])/length(X3))^2
    Xt3<-(length(X3[X3==2])/length(X3))^2
    Xg3<-(length(X3[X3==3])/length(X3))^2
    Xc3<-(length(X3[X3==4])/length(X3))^2
    Xx3<-(length(X3[X3==5])/length(X3))^2
    SX3<-1-sum(Xa3,Xt3,Xg3,Xc3,Xx3)
    V4<-c(Xa3,Xt3,Xg3,Xc3,Xx3)
    
    X4<-as.vector(Depth2[rownames(CO),i])
    Xa4<-(length(X4[X4==1])/length(X4))^2
    Xt4<-(length(X4[X4==2])/length(X4))^2
    Xg4<-(length(X4[X4==3])/length(X4))^2
    Xc4<-(length(X4[X4==4])/length(X4))^2
    Xx4<-(length(X4[X4==5])/length(X4))^2
    SX4<-1-sum(Xa4,Xt4,Xg4,Xc4,Xx4)
    V5<-c(Xa4,Xt4,Xg4,Xc4,Xx4)
    
    C12id<-cor(V1,V2)
    C13id<-cor(V1,V3)
    C14id<-cor(V1,V4)
    C15id<-cor(V1,V5)
    
    C21mt<-cor(V2,V1)
    C23mt<-cor(V2,V3)
    C24mt<-cor(V2,V4)
    C25mt<-cor(V2,V5)
    
    C31or<-cor(V3,V1)
    C32or<-cor(V3,V2)
    C34or<-cor(V3,V4)
    C35or<-cor(V3,V5)
    
    C41wy<-cor(V4,V1)
    C42wy<-cor(V4,V2)
    C43wy<-cor(V4,V3)
    C45wy<-cor(V4,V5)
    
    C51co<-cor(V5,V1)
    C52co<-cor(V5,V2)
    C53co<-cor(V5,V3)
    C54co<-cor(V5,V4)
    
    C1cor<-c(C12id,C13id,C14id,C15id,C21mt,C23mt,C24mt,C25mt,
             C31or,C32or,C34or,C35or,C41wy,C42wy,C43wy,C45wy,
             C51co,C52co,C53co,C54co)
    
    AFV1<-c(SX,SX1,SX2,SX3,SX4)
    AFV<-rbind(AFV,AFV1)
    Cor1<-rbind(Cor1,C1cor)
  }
  AFV<-AFV[-1,]
  Cor1<-(1-abs(Cor1))
  
  
  Cor1
  COcor<-rowMeans(Cor1[,c(1,2,3,4)])
  IDcor<-rowMeans(Cor1[,c(5,6,7,8)])
  NEcor<-rowMeans(Cor1[,c(9,10,11,12)])
  ORcor<-rowMeans(Cor1[,c(12,13,14,15)])
  CAcor<-rowMeans(Cor1[,c(16,17,18,19)])
  
  ORcor[ORcor<0.5]=0
  IDcor[IDcor<0.5]=0
  NEcor[NEcor<0.5]=0
  ORcor[ORcor<0.5]=0
  CAcor[COcor<0.5]=0
  
  par(mfrow=c(1,1))
  plot(NULL,ylim=c(0,40),xlim = c(0,3175))
  
  
  ORFs<-read.table('out.cds.blast')
  ORFF<-ORFs[ORFs[,1]=='giCO',c(13,14)]
  
  lorf<-c()
  for (i in seq(1,length(ORFF[,1]))){
    print(max(c(ORFF[i,1],ORFF[i,2]))-min(c(ORFF[i,1],ORFF[i,2])))
    lorf<-c(lorf,max(c(ORFF[i,1],ORFF[i,2]))-min(c(ORFF[i,1],ORFF[i,2])))
  }
  olorf<-order(lorf)
  oolorf<-olorf[seq(1,length(lorf)-7)]
  ORFF<-ORFF[oolorf,]
  
  
  lorf<-c()
  for (p in seq(1,length(ORFF[,1]))){
    ll<-length(seq(ORFF[p,1],ORFF[p,2]))
    lorf<-c(lorf,ll)
  }
  
  lorf<-c()
  for (p in seq(1,length(ORFF[,1]))){
    ll<-(length(seq(ORFF[p,1],ORFF[p,2])))
    lorf<-c(lorf,ll)
    
  }
  lorf
  lorfo<-order(lorf,decreasing = F)
  
  Pdf<-c(rep(0,3175))
  Pfd1<-c(0,0,0)
  for (p in lorfo){
    print(ORFF[p,])
    P1<-c(rep(0,3175))
    P1[seq(ORFF[p,1],ORFF[p,2])]=1
    Pdf<-rbind(Pdf,P1)
    print((colSums(Pdf[,seq(ORFF[p,1],ORFF[p,2])])))
    MCS<-max(colSums(Pdf[,seq(ORFF[p,1],ORFF[p,2])]))
    print(MCS)
    ORFFF<-(cbind(ORFF[p,],MCS))
    Pfd1<-rbind(Pfd1,ORFFF)
  }
  
  Pfd1<-Pfd1[-1,]
  inens<-c()
  for (p in seq(1,length(Pfd1[,1]))){
    inens<-c(inens,1)
  }
  
  
  for (p in seq(1,length(Pfd1[,1]))){
    color_transparent <- adjustcolor('grey', alpha.f = inens[p])
    xbottom <- Pfd1[p,1]
    ybottom <- Pfd1[p,3]+19
    xtop <- Pfd1[p,2]
    ytop<- Pfd1[p,3]+20
    xpol <- c(xbottom,xbottom,xtop,xtop)
    ypol <- c(ybottom,ytop,ytop,ybottom)
    polygon(x=cuts(t(xpol)), y=cuts(t(ypol)), col = color_transparent, border ='black')
  }
  
  #ADD sncRNA
  sncRNA<-read.table('out.snc.blast')
  Pfd1<-sncRNA[,c(9,10)]
  for (p in seq(1,length(Pfd1[,1]))){
    color_transparent <- adjustcolor('black', alpha.f = .5)
    xbottom <- Pfd1[p,1]
    ybottom <- 18
    xtop <- Pfd1[p,2]
    ytop<- 19
    xpol <- c(xbottom,xbottom,xtop,xtop)
    ypol <- c(ybottom,ytop,ytop,ybottom)
    polygon(x=cuts(t(xpol)), y=cuts(t(ypol)), col = color_transparent, border =NA)
  }
  
  
  
  for (p in seq(1,length(ORcor))){
    inens<- ORcor[p]*2
    color_transparent <- adjustcolor('orange', alpha.f = inens)
    Gm<-p
    GM<-p+1
    xbottom <- Gm
    ybottom <- 17
    xtop <- GM
    ytop<- 18
    xpol <- c(xbottom,xbottom,xtop,xtop)
    ypol <- c(ybottom,ytop,ytop,ybottom)
    polygon(x=cuts(t(xpol)), y=cuts(t(ypol)), col = color_transparent, border =NA)
  }
  
  for (p in seq(1,length(IDcor))){
    inens<- IDcor[p]*2
    color_transparent <- adjustcolor('blue', alpha.f = inens)
    Gm<-p
    GM<-p+1
    xbottom <- Gm
    ybottom <- 16
    xtop <- GM
    ytop<- 17
    xpol <- c(xbottom,xbottom,xtop,xtop)
    ypol <- c(ybottom,ytop,ytop,ybottom)
    polygon(x=cuts(t(xpol)), y=cuts(t(ypol)), col = color_transparent, border =NA)
  }
  
  for (p in seq(1,length(NEcor))){
    inens<- NEcor[p]*2
    color_transparent <- adjustcolor('green', alpha.f = inens)
    Gm<-p
    GM<-p+1
    xbottom <- Gm
    ybottom <- 15
    xtop <- GM
    ytop<- 16
    xpol <- c(xbottom,xbottom,xtop,xtop)
    ypol <- c(ybottom,ytop,ytop,ybottom)
    polygon(x=cuts(t(xpol)), y=cuts(t(ypol)), col = color_transparent, border =NA)
  }
  
  for (p in seq(1,length(CAcor))){
    inens<- CAcor[p]*2
    color_transparent <- adjustcolor('red', alpha.f = inens)
    Gm<-p
    GM<-p+1
    xbottom <- Gm
    ybottom <- 14
    xtop <- GM
    ytop<- 15
    xpol <- c(xbottom,xbottom,xtop,xtop)
    ypol <- c(ybottom,ytop,ytop,ybottom)
    polygon(x=cuts(t(xpol)), y=cuts(t(ypol)), col = color_transparent, border =NA)
  }
  
  for (p in seq(1,length(COcor))){
    inens<- COcor[p]*2
    color_transparent <- adjustcolor('purple', alpha.f = inens)
    Gm<-p
    GM<-p+1
    xbottom <- Gm
    ybottom <- 13
    xtop <- GM
    ytop<- 14
    xpol <- c(xbottom,xbottom,xtop,xtop)
    ypol <- c(ybottom,ytop,ytop,ybottom)
    polygon(x=cuts(t(xpol)), y=cuts(t(ypol)), col = color_transparent, border =NA)
  }
  
  
  cbind(ORcor,IDcor,NEcor,CAcor,COcor)
  
  
  
  Depth2<-Depth1[rownames(SEVRN),]
  rownames(Depth2)
  SEVRN[,4]
  Depth2<-matrix(as.numeric(Depth2),ncol = 3175)
  
  d <- dist(Depth2) # euclidean distances between the rows
  fit <- cmdscale(d,eig=TRUE, k=2) # k is the number of dim
  fit # view results
  
  # plot solution
  x <- fit$points[,1]
  y <- fit$points[,2]
  
  SI<-as.numeric(SEVRN[,4])
  SI
  unique(SI)
  orderN<-order(SI)
  color_colors<-rainbow(13)
  
  plot(x, y, xlab="Coordinate 1", ylab="Coordinate 2",
       main="BCTV Virus Genome Sequences", type="p",
       pch=19,col=color_colors[SI])
  
  legend(x = "topright",       # Position
         legend = unique(SEVRN[,4]),  # Legend texts
         pch = 19,           # Line types
         col = color_colors[unique(SI)],           # Line colors
         lwd = 0)                 # Line width
  
  
  
  library(ade4)
  unique(DepthMD[,2])
  unique(DepthMD[,2])
  ROWNamova<-rownames(DepthMD[DepthMD[,2]=='BCTV-CO',])
  Structure<-as.data.frame(as.vector(DepthMD[DepthMD[,2]=='BCTV-CO',4]))
  colnames(Structure)<-'State'
  SAMP<-t(as.data.frame(Depth2[ROWNamova,]))
  colnames(SAMP)<-seq(1,length(SAMP[1,]))
  SAMP1<-cbind(SAMP,AFV[,4])
  SAMP2<-SAMP1[SAMP1[,20]>0,]
  SAMP3<-SAMP2[,-112]
  SAMP4<-as.matrix(SAMP3)
  SAMP5<-c()
  for (i in seq(1,length(SAMP4[,1]))){
    if (length(unique(SAMP4[i,]))>1){
      SAMP5<-rbind(SAMP5,SAMP4[i,])
    }
  }
  
  mean(rowMeans(Cor1[,c(1,2,3)]))
  mean(rowMeans(Cor1[,c(4,5,6)]))
  mean(rowMeans(Cor1[,c(7,8,9)]))
  mean(rowMeans(Cor1[,c(10,11,12)]))
  hist(rowMeans(Cor1[,c(1,2,3)]))
  hist(rowMeans(Cor1[,c(4,5,6)]))
  hist(rowMeans(Cor1[,c(7,8,9)]))
  hist(rowMeans(Cor1[,c(10,11,12)]))
  
  VectorRM<-rowMeans(Cor1[,c(1,2,3)])
  VectorRM<-data.frame(VectorRM)
  VectorRM<-cbind(VectorRM,VectorRM)
  VectorRM[VectorRM[,1]>0.5,]
  length(VectorRM[VectorRM[,1]>0.5,1])
  
  
  AMOV<-amova(samples=data.frame(SAMP5),structures = Structure)
  AMOV$results
  AMOV$componentsofcovariance
  AMOV$statphi
  
  for (i in seq(1,length(SAMP5[,1]))){
    AMOV<-amova(samples = data.frame(t(SAMP5[i,])),structures = Structure)
    print(AMOV$componentsofcovariance[1,])
  }
  
  library(dendextend)
  hc<-hclust(DIST)
  dhc <- as.dendrogram(hc)
  plot(dhc)
  dend2 <- cut(dhc, h = 17)
  dend2 <- color_branches(dhc, k = 12)
  dend2<-set(dend2, "branches_lwd", 3)
  plot(dend2)
  
  # creates a single item vector of the clusters    
  myclusters <- cutree(dhc, k=12)
  
  # make the dataframe of two columns cluster number and label
  clusterDF <-  data.frame(Cluster = as.numeric(unlist(myclusters)),
                           Branch = names(myclusters))
  clusterDF[order(clusterDF[,1]),]
  ROWNamova[clusterDF[order(clusterDF[,1]),2]]
  
  # sort by cluster ascending
  clusterDFSort <- clusterDF %>% arrange(Cluster)
  Structure
  
}

unique(WORRN[,4]) 
#WORLAND
{
  unique(WORRN[,4])
  
  FstALL<-data.frame()
  Cor1<-data.frame()
  AFV<-data.frame(1,2,3,4,5,6)
  colnames(AFV)<-c('IR','MT','OR','ID','WY','CO')
  ID<-SEVRN[SEVRN[,4]=='ID',]
  MT<-SEVRN[SEVRN[,4]=='MT',]
  OR<-SEVRN[SEVRN[,4]=='OR',]
  WY<-SEVRN[SEVRN[,4]=='WY',]
  CO<-SEVRN[SEVRN[,4]=='CO',]
  
  for (i in seq(1,3175)){
    i<-i
    print(i)
    X<-as.vector(Depth2[rownames(ID),i])
    Xa<-(length(X[X==1])/length(X))^2
    Xt<-(length(X[X==2])/length(X))^2
    Xg<-(length(X[X==3])/length(X))^2
    Xc<-(length(X[X==4])/length(X))^2
    Xx<-(length(X[X==5])/length(X))^2
    SX<-1-sum(Xa,Xt,Xg,Xc,Xx)
    V1<-c(Xa,Xt,Xg,Xc,Xx)
    SUMX<-(SX)
    
    X1<-as.vector(Depth2[rownames(MT),i])
    Xa1<-(length(X1[X1==1])/length(X1))^2
    Xt1<-(length(X1[X1==2])/length(X1))^2
    Xg1<-(length(X1[X1==3])/length(X1))^2
    Xc1<-(length(X1[X1==4])/length(X1))^2
    Xx1<-(length(X1[X1==5])/length(X1))^2
    SX1<-1-sum(Xa1,Xt1,Xg1,Xc1,Xx1)
    SUMX1<-(SX1)
    FstSite<-((SX+SX1)-(SX1))/(SX+SX1)
    V2<-c(Xa1,Xt1,Xg1,Xc1,Xx1)
    C1<-cor(V1,V2)
    
    X2<-as.vector(Depth2[rownames(OR),i])
    Xa2<-(length(X2[X2==1])/length(X2))^2
    Xt2<-(length(X2[X2==2])/length(X2))^2
    Xg2<-(length(X2[X2==3])/length(X2))^2
    Xc2<-(length(X2[X2==4])/length(X2))^2
    Xx2<-(length(X2[X2==5])/length(X2))^2
    SX2<-1-sum(Xa2,Xt2,Xg2,Xc2,Xx2)
    V3<-c(Xa2,Xt2,Xg2,Xc2,Xx2)
    
    X3<-as.vector(Depth2[rownames(WY),i])
    Xa3<-(length(X3[X3==1])/length(X3))^2
    Xt3<-(length(X3[X3==2])/length(X3))^2
    Xg3<-(length(X3[X3==3])/length(X3))^2
    Xc3<-(length(X3[X3==4])/length(X3))^2
    Xx3<-(length(X3[X3==5])/length(X3))^2
    SX3<-1-sum(Xa3,Xt3,Xg3,Xc3,Xx3)
    V4<-c(Xa3,Xt3,Xg3,Xc3,Xx3)
    
    X4<-as.vector(Depth2[rownames(CO),i])
    Xa4<-(length(X4[X4==1])/length(X4))^2
    Xt4<-(length(X4[X4==2])/length(X4))^2
    Xg4<-(length(X4[X4==3])/length(X4))^2
    Xc4<-(length(X4[X4==4])/length(X4))^2
    Xx4<-(length(X4[X4==5])/length(X4))^2
    SX4<-1-sum(Xa4,Xt4,Xg4,Xc4,Xx4)
    V5<-c(Xa4,Xt4,Xg4,Xc4,Xx4)
    
    C12id<-cor(V1,V2)
    C13id<-cor(V1,V3)
    C14id<-cor(V1,V4)
    C15id<-cor(V1,V5)
    
    C21mt<-cor(V2,V1)
    C23mt<-cor(V2,V3)
    C24mt<-cor(V2,V4)
    C25mt<-cor(V2,V5)
    
    C31or<-cor(V3,V1)
    C32or<-cor(V3,V2)
    C34or<-cor(V3,V4)
    C35or<-cor(V3,V5)
    
    C41wy<-cor(V4,V1)
    C42wy<-cor(V4,V2)
    C43wy<-cor(V4,V3)
    C45wy<-cor(V4,V5)
    
    C51co<-cor(V5,V1)
    C52co<-cor(V5,V2)
    C53co<-cor(V5,V3)
    C54co<-cor(V5,V4)
    
    C1cor<-c(C12id,C13id,C14id,C15id,C21mt,C23mt,C24mt,C25mt,
             C31or,C32or,C34or,C35or,C41wy,C42wy,C43wy,C45wy,
             C51co,C52co,C53co,C54co)
    
    AFV1<-c(SX,SX1,SX2,SX3,SX4)
    AFV<-rbind(AFV,AFV1)
    Cor1<-rbind(Cor1,C1cor)
  }
  AFV<-AFV[-1,]
  Cor1<-(1-abs(Cor1))
  
  
  Cor1
  IDcor<-rowMeans(Cor1[,c(1,2,3,4)])
  MTcor<-rowMeans(Cor1[,c(5,6,7,8)])
  ORcor<-rowMeans(Cor1[,c(9,10,11,12)])
  WYcor<-rowMeans(Cor1[,c(12,13,14,15)])
  COcor<-rowMeans(Cor1[,c(16,17,18,19)])
  
  ORcor[ORcor<0.5]=0
  IDcor[IDcor<0.5]=0
  MTcor[MTcor<0.5]=0
  WYcor[WYcor<0.5]=0
  COcor[COcor<0.5]=0
  
  par(mfrow=c(1,1))
  plot(NULL,ylim=c(0,40),xlim = c(0,3175))
  
  ORFF<-ORFs[ORFs[,1]=='giSEV',c(13,14)]
  
  lorf<-c()
  for (i in seq(1,length(ORFF[,1]))){
    print(max(c(ORFF[i,1],ORFF[i,2]))-min(c(ORFF[i,1],ORFF[i,2])))
    lorf<-c(lorf,max(c(ORFF[i,1],ORFF[i,2]))-min(c(ORFF[i,1],ORFF[i,2])))
  }
  olorf<-order(lorf)
  oolorf<-olorf[seq(1,length(lorf)-7)]
  ORFF<-ORFF[oolorf,]
  
  
  lorf<-c()
  for (p in seq(1,length(ORFF[,1]))){
    ll<-length(seq(ORFF[p,1],ORFF[p,2]))
    lorf<-c(lorf,ll)
  }
  
  lorf<-c()
  for (p in seq(1,length(ORFF[,1]))){
    ll<-(length(seq(ORFF[p,1],ORFF[p,2])))
    lorf<-c(lorf,ll)
    
  }
  lorf
  lorfo<-order(lorf,decreasing = F)
  
  Pdf<-c(rep(0,3175))
  Pfd1<-c(0,0,0)
  for (p in lorfo){
    print(ORFF[p,])
    P1<-c(rep(0,3175))
    P1[seq(ORFF[p,1],ORFF[p,2])]=1
    Pdf<-rbind(Pdf,P1)
    print((colSums(Pdf[,seq(ORFF[p,1],ORFF[p,2])])))
    MCS<-max(colSums(Pdf[,seq(ORFF[p,1],ORFF[p,2])]))
    print(MCS)
    ORFFF<-(cbind(ORFF[p,],MCS))
    Pfd1<-rbind(Pfd1,ORFFF)
  }
  
  Pfd1<-Pfd1[-1,]
  inens<-c()
  for (p in seq(1,length(Pfd1[,1]))){
    inens<-c(inens,1)
  }
  
  
  for (p in seq(1,length(Pfd1[,1]))){
    color_transparent <- adjustcolor('grey', alpha.f = inens[p])
    xbottom <- Pfd1[p,1]
    ybottom <- Pfd1[p,3]+19
    xtop <- Pfd1[p,2]
    ytop<- Pfd1[p,3]+20
    xpol <- c(xbottom,xbottom,xtop,xtop)
    ypol <- c(ybottom,ytop,ytop,ybottom)
    polygon(x=cuts(t(xpol)), y=cuts(t(ypol)), col = color_transparent, border ='black')
  }
  
  
  #ADD sncRNA
  sncRNA<-read.table('out.snc.blast')
  Pfd1<-sncRNA[,c(9,10)]
  for (p in seq(1,length(Pfd1[,1]))){
    color_transparent <- adjustcolor('black', alpha.f = .5)
    xbottom <- Pfd1[p,1]
    ybottom <- 18
    xtop <- Pfd1[p,2]
    ytop<- 19
    xpol <- c(xbottom,xbottom,xtop,xtop)
    ypol <- c(ybottom,ytop,ytop,ybottom)
    polygon(x=cuts(t(xpol)), y=cuts(t(ypol)), col = color_transparent, border =NA)
  }
  
  
  
  for (p in seq(1,length(ORcor))){
    inens<- ORcor[p]*2
    color_transparent <- adjustcolor('orange', alpha.f = inens)
    Gm<-p
    GM<-p+1
    xbottom <- Gm
    ybottom <- 17
    xtop <- GM
    ytop<- 18
    xpol <- c(xbottom,xbottom,xtop,xtop)
    ypol <- c(ybottom,ytop,ytop,ybottom)
    polygon(x=cuts(t(xpol)), y=cuts(t(ypol)), col = color_transparent, border =NA)
  }
  
  for (p in seq(1,length(IDcor))){
    inens<- IDcor[p]*2
    color_transparent <- adjustcolor('blue', alpha.f = inens)
    Gm<-p
    GM<-p+1
    xbottom <- Gm
    ybottom <- 16
    xtop <- GM
    ytop<- 17
    xpol <- c(xbottom,xbottom,xtop,xtop)
    ypol <- c(ybottom,ytop,ytop,ybottom)
    polygon(x=cuts(t(xpol)), y=cuts(t(ypol)), col = color_transparent, border =NA)
  }
  
  for (p in seq(1,length(MTcor))){
    inens<- MTcor[p]*2
    color_transparent <- adjustcolor('green', alpha.f = inens)
    Gm<-p
    GM<-p+1
    xbottom <- Gm
    ybottom <- 15
    xtop <- GM
    ytop<- 16
    xpol <- c(xbottom,xbottom,xtop,xtop)
    ypol <- c(ybottom,ytop,ytop,ybottom)
    polygon(x=cuts(t(xpol)), y=cuts(t(ypol)), col = color_transparent, border =NA)
  }
  
  for (p in seq(1,length(WYcor))){
    inens<- WYcor[p]*2
    color_transparent <- adjustcolor('red', alpha.f = inens)
    Gm<-p
    GM<-p+1
    xbottom <- Gm
    ybottom <- 14
    xtop <- GM
    ytop<- 15
    xpol <- c(xbottom,xbottom,xtop,xtop)
    ypol <- c(ybottom,ytop,ytop,ybottom)
    polygon(x=cuts(t(xpol)), y=cuts(t(ypol)), col = color_transparent, border =NA)
  }
  
  for (p in seq(1,length(COcor))){
    inens<- COcor[p]*2
    color_transparent <- adjustcolor('purple', alpha.f = inens)
    Gm<-p
    GM<-p+1
    xbottom <- Gm
    ybottom <- 13
    xtop <- GM
    ytop<- 14
    xpol <- c(xbottom,xbottom,xtop,xtop)
    ypol <- c(ybottom,ytop,ytop,ybottom)
    polygon(x=cuts(t(xpol)), y=cuts(t(ypol)), col = color_transparent, border =NA)
  }
  
  
  cbind(ORcor,IDcor,MTcor,WYcor,COcor)
  
  
  
  Depth2<-Depth1[rownames(SEVRN),]
  rownames(Depth2)
  SEVRN[,4]
  Depth2<-matrix(as.numeric(Depth2),ncol = 3175)
  
  d <- dist(Depth2) # euclidean distances between the rows
  fit <- cmdscale(d,eig=TRUE, k=2) # k is the number of dim
  fit # view results
  
  # plot solution
  x <- fit$points[,1]
  y <- fit$points[,2]
  
  SI<-as.numeric(SEVRN[,4])
  SI
  unique(SI)
  orderN<-order(SI)
  color_colors<-rainbow(13)
  
  plot(x, y, xlab="Coordinate 1", ylab="Coordinate 2",
       main="BCTV Virus Genome Sequences", type="p",
       pch=19,col=color_colors[SI])
  
  legend(x = "topright",       # Position
         legend = unique(SEVRN[,4]),  # Legend texts
         pch = 19,           # Line types
         col = color_colors[unique(SI)],           # Line colors
         lwd = 0)                 # Line width
  
}

#SEVERE
{
  unique(SEVRN[,4])
  unique(DepthMD[,2])
  ROWNamova<-rownames(DepthMD[DepthMD[,2]=='BCTV-Svr',])
  Structure<-as.data.frame(as.vector(DepthMD[DepthMD[,2]=='BCTV-Svr',4]))
  colnames(Structure)<-'State'
  SAMP<-t(as.data.frame(Depth2[ROWNamova,]))
  colnames(SAMP)<-seq(1,length(SAMP[1,]))
  SAMP1<-cbind(SAMP,AFV[,4])
  SAMP2<-SAMP1[SAMP1[,20]>0,]
  SAMP3<-SAMP2[,-25]
  SAMP4<-as.matrix(SAMP3)
  SAMP5<-c()
  for (i in seq(1,length(SAMP4[,1]))){
    if (length(unique(SAMP4[i,]))>1){
      SAMP5<-rbind(SAMP5,SAMP4[i,])
    }
  }
  
  SMAP5<-data.frame(SAMP5)
  DIST<-dist(t(SAMP5),method='euclidean')
  library(dendextend)
  hc<-hclust(DIST)
  dhc <- as.dendrogram(hc)
  dend2 <- color_branches(dhc, k = 6)
  dend2<-set(dend2, "branches_lwd", 3)
  plot(dend2)
  # creates a single item vector of the clusters    
  myclusters <- cutree(dhc, k=6)
  # make the dataframe of two columns cluster number and label
  clusterDF <-  data.frame(Cluster = as.numeric(unlist(myclusters)),
                           Branch = names(myclusters))
  clusterDF[order(clusterDF[,1]),]
  ROWNamova[clusterDF[order(clusterDF[,1]),2]]
  cbind(DepthMD[ROWNamova[clusterDF[order(clusterDF[,1]),2]],],clusterDF[order(clusterDF[,1]),])
  
  
  
  FstALL<-data.frame()
  Cor1<-data.frame()
  AFV<-data.frame(1,2,3,4,5,6)
  colnames(AFV)<-c('IR','MT','OR','ID','WY','CO')
  ID<-SEVRN[SEVRN[,4]=='ID',]
  MT<-SEVRN[SEVRN[,4]=='MT',]
  OR<-SEVRN[SEVRN[,4]=='OR',]
  WY<-SEVRN[SEVRN[,4]=='WY',]
  CO<-SEVRN[SEVRN[,4]=='CO',]
  
  for (i in seq(1,3175)){
    i<-i
    print(i)
    X<-as.vector(Depth2[rownames(ID),i])
    Xa<-(length(X[X==1])/length(X))^2
    Xt<-(length(X[X==2])/length(X))^2
    Xg<-(length(X[X==3])/length(X))^2
    Xc<-(length(X[X==4])/length(X))^2
    Xx<-(length(X[X==5])/length(X))^2
    SX<-1-sum(Xa,Xt,Xg,Xc,Xx)
    V1<-c(Xa,Xt,Xg,Xc,Xx)
    SUMX<-(SX)
    
    X1<-as.vector(Depth2[rownames(MT),i])
    Xa1<-(length(X1[X1==1])/length(X1))^2
    Xt1<-(length(X1[X1==2])/length(X1))^2
    Xg1<-(length(X1[X1==3])/length(X1))^2
    Xc1<-(length(X1[X1==4])/length(X1))^2
    Xx1<-(length(X1[X1==5])/length(X1))^2
    SX1<-1-sum(Xa1,Xt1,Xg1,Xc1,Xx1)
    SUMX1<-(SX1)
    FstSite<-((SX+SX1)-(SX1))/(SX+SX1)
    V2<-c(Xa1,Xt1,Xg1,Xc1,Xx1)
    C1<-cor(V1,V2)
    
    X2<-as.vector(Depth2[rownames(OR),i])
    Xa2<-(length(X2[X2==1])/length(X2))^2
    Xt2<-(length(X2[X2==2])/length(X2))^2
    Xg2<-(length(X2[X2==3])/length(X2))^2
    Xc2<-(length(X2[X2==4])/length(X2))^2
    Xx2<-(length(X2[X2==5])/length(X2))^2
    SX2<-1-sum(Xa2,Xt2,Xg2,Xc2,Xx2)
    V3<-c(Xa2,Xt2,Xg2,Xc2,Xx2)
    
    X3<-as.vector(Depth2[rownames(WY),i])
    Xa3<-(length(X3[X3==1])/length(X3))^2
    Xt3<-(length(X3[X3==2])/length(X3))^2
    Xg3<-(length(X3[X3==3])/length(X3))^2
    Xc3<-(length(X3[X3==4])/length(X3))^2
    Xx3<-(length(X3[X3==5])/length(X3))^2
    SX3<-1-sum(Xa3,Xt3,Xg3,Xc3,Xx3)
    V4<-c(Xa3,Xt3,Xg3,Xc3,Xx3)
    
    X4<-as.vector(Depth2[rownames(CO),i])
    Xa4<-(length(X4[X4==1])/length(X4))^2
    Xt4<-(length(X4[X4==2])/length(X4))^2
    Xg4<-(length(X4[X4==3])/length(X4))^2
    Xc4<-(length(X4[X4==4])/length(X4))^2
    Xx4<-(length(X4[X4==5])/length(X4))^2
    SX4<-1-sum(Xa4,Xt4,Xg4,Xc4,Xx4)
    V5<-c(Xa4,Xt4,Xg4,Xc4,Xx4)
    
    C12id<-cor(V1,V2)
    C13id<-cor(V1,V3)
    C14id<-cor(V1,V4)
    C15id<-cor(V1,V5)
    
    C21mt<-cor(V2,V1)
    C23mt<-cor(V2,V3)
    C24mt<-cor(V2,V4)
    C25mt<-cor(V2,V5)
    
    C31or<-cor(V3,V1)
    C32or<-cor(V3,V2)
    C34or<-cor(V3,V4)
    C35or<-cor(V3,V5)
    
    C41wy<-cor(V4,V1)
    C42wy<-cor(V4,V2)
    C43wy<-cor(V4,V3)
    C45wy<-cor(V4,V5)
    
    C51co<-cor(V5,V1)
    C52co<-cor(V5,V2)
    C53co<-cor(V5,V3)
    C54co<-cor(V5,V4)
    
    C1cor<-c(C12id,C13id,C14id,C15id,C21mt,C23mt,C24mt,C25mt,
             C31or,C32or,C34or,C35or,C41wy,C42wy,C43wy,C45wy,
             C51co,C52co,C53co,C54co)
    
    AFV1<-c(SX,SX1,SX2,SX3,SX4)
    AFV<-rbind(AFV,AFV1)
    Cor1<-rbind(Cor1,C1cor)
  }
  AFV<-AFV[-1,]
  Cor1<-(1-abs(Cor1))
  
  
  Cor1
  IDcor<-rowMeans(Cor1[,c(1,2,3,4)])
  MTcor<-rowMeans(Cor1[,c(5,6,7,8)])
  ORcor<-rowMeans(Cor1[,c(9,10,11,12)])
  WYcor<-rowMeans(Cor1[,c(12,13,14,15)])
  COcor<-rowMeans(Cor1[,c(16,17,18,19)])
  
  ORcor[ORcor<0.5]=0
  IDcor[IDcor<0.5]=0
  MTcor[MTcor<0.5]=0
  WYcor[WYcor<0.5]=0
  COcor[COcor<0.5]=0
  
  par(mfrow=c(1,1))
  plot(NULL,ylim=c(0,40),xlim = c(0,3175))
  
  ORFF<-ORFs[ORFs[,1]=='giSEV',c(13,14)]
  
  lorf<-c()
  for (i in seq(1,length(ORFF[,1]))){
    print(max(c(ORFF[i,1],ORFF[i,2]))-min(c(ORFF[i,1],ORFF[i,2])))
    lorf<-c(lorf,max(c(ORFF[i,1],ORFF[i,2]))-min(c(ORFF[i,1],ORFF[i,2])))
  }
  olorf<-order(lorf)
  oolorf<-olorf[seq(1,length(lorf)-7)]
  ORFF<-ORFF[oolorf,]
  
  
  lorf<-c()
  for (p in seq(1,length(ORFF[,1]))){
    ll<-length(seq(ORFF[p,1],ORFF[p,2]))
    lorf<-c(lorf,ll)
  }
  
  lorf<-c()
  for (p in seq(1,length(ORFF[,1]))){
    ll<-(length(seq(ORFF[p,1],ORFF[p,2])))
    lorf<-c(lorf,ll)
    
  }
  lorf
  lorfo<-order(lorf,decreasing = F)
  
  Pdf<-c(rep(0,3175))
  Pfd1<-c(0,0,0)
  for (p in lorfo){
    print(ORFF[p,])
    P1<-c(rep(0,3175))
    P1[seq(ORFF[p,1],ORFF[p,2])]=1
    Pdf<-rbind(Pdf,P1)
    print((colSums(Pdf[,seq(ORFF[p,1],ORFF[p,2])])))
    MCS<-max(colSums(Pdf[,seq(ORFF[p,1],ORFF[p,2])]))
    print(MCS)
    ORFFF<-(cbind(ORFF[p,],MCS))
    Pfd1<-rbind(Pfd1,ORFFF)
  }
  
  Pfd1<-Pfd1[-1,]
  inens<-c()
  for (p in seq(1,length(Pfd1[,1]))){
    inens<-c(inens,1)
  }
  
  
  for (p in seq(1,length(Pfd1[,1]))){
    color_transparent <- adjustcolor('grey', alpha.f = inens[p])
    xbottom <- Pfd1[p,1]
    ybottom <- Pfd1[p,3]+19
    xtop <- Pfd1[p,2]
    ytop<- Pfd1[p,3]+20
    xpol <- c(xbottom,xbottom,xtop,xtop)
    ypol <- c(ybottom,ytop,ytop,ybottom)
    polygon(x=cuts(t(xpol)), y=cuts(t(ypol)), col = color_transparent, border ='black')
  }
  
  
  #ADD sncRNA
  sncRNA<-read.table('out.snc.blast')
  Pfd1<-sncRNA[,c(9,10)]
  for (p in seq(1,length(Pfd1[,1]))){
    color_transparent <- adjustcolor('black', alpha.f = .5)
    xbottom <- Pfd1[p,1]
    ybottom <- 18
    xtop <- Pfd1[p,2]
    ytop<- 19
    xpol <- c(xbottom,xbottom,xtop,xtop)
    ypol <- c(ybottom,ytop,ytop,ybottom)
    polygon(x=cuts(t(xpol)), y=cuts(t(ypol)), col = color_transparent, border =NA)
  }
  
  
  
  for (p in seq(1,length(ORcor))){
    inens<- ORcor[p]*2
    color_transparent <- adjustcolor('orange', alpha.f = inens)
    Gm<-p
    GM<-p+1
    xbottom <- Gm
    ybottom <- 17
    xtop <- GM
    ytop<- 18
    xpol <- c(xbottom,xbottom,xtop,xtop)
    ypol <- c(ybottom,ytop,ytop,ybottom)
    polygon(x=cuts(t(xpol)), y=cuts(t(ypol)), col = color_transparent, border =NA)
  }
  
  for (p in seq(1,length(IDcor))){
    inens<- IDcor[p]*2
    color_transparent <- adjustcolor('blue', alpha.f = inens)
    Gm<-p
    GM<-p+1
    xbottom <- Gm
    ybottom <- 16
    xtop <- GM
    ytop<- 17
    xpol <- c(xbottom,xbottom,xtop,xtop)
    ypol <- c(ybottom,ytop,ytop,ybottom)
    polygon(x=cuts(t(xpol)), y=cuts(t(ypol)), col = color_transparent, border =NA)
  }
  
  for (p in seq(1,length(MTcor))){
    inens<- MTcor[p]*2
    color_transparent <- adjustcolor('green', alpha.f = inens)
    Gm<-p
    GM<-p+1
    xbottom <- Gm
    ybottom <- 15
    xtop <- GM
    ytop<- 16
    xpol <- c(xbottom,xbottom,xtop,xtop)
    ypol <- c(ybottom,ytop,ytop,ybottom)
    polygon(x=cuts(t(xpol)), y=cuts(t(ypol)), col = color_transparent, border =NA)
  }
  
  for (p in seq(1,length(WYcor))){
    inens<- WYcor[p]*2
    color_transparent <- adjustcolor('red', alpha.f = inens)
    Gm<-p
    GM<-p+1
    xbottom <- Gm
    ybottom <- 14
    xtop <- GM
    ytop<- 15
    xpol <- c(xbottom,xbottom,xtop,xtop)
    ypol <- c(ybottom,ytop,ytop,ybottom)
    polygon(x=cuts(t(xpol)), y=cuts(t(ypol)), col = color_transparent, border =NA)
  }
  
  for (p in seq(1,length(COcor))){
    inens<- COcor[p]*2
    color_transparent <- adjustcolor('purple', alpha.f = inens)
    Gm<-p
    GM<-p+1
    xbottom <- Gm
    ybottom <- 13
    xtop <- GM
    ytop<- 14
    xpol <- c(xbottom,xbottom,xtop,xtop)
    ypol <- c(ybottom,ytop,ytop,ybottom)
    polygon(x=cuts(t(xpol)), y=cuts(t(ypol)), col = color_transparent, border =NA)
  }
  
}
  
  
  
  
  cbind(ORcor,IDcor,MTcor,WYcor,COcor)
  
  Depth2<-Depth1[rownames(SEVRN),]
  rownames(Depth2)
  SEVRN[,4]
  Depth2<-matrix(as.numeric(Depth2),ncol = 3175)
  
  d <- dist(Depth2) # euclidean distances between the rows
  fit <- cmdscale(d,eig=TRUE, k=2) # k is the number of dim
  fit # view results
  
  # plot solution
  x <- fit$points[,1]
  y <- fit$points[,2]
  
  SI<-as.numeric(SEVRN[,4])
  SI
  unique(SI)
  orderN<-order(SI)
  color_colors<-rainbow(13)
  
  plot(x, y, xlab="Coordinate 1", ylab="Coordinate 2",
       main="BCTV Virus Genome Sequences", type="p",
       pch=19,col=color_colors[SI])
  
  legend(x = "topright",       # Position
         legend = unique(SEVRN[,4]),  # Legend texts
         pch = 19,           # Line types
         col = color_colors[unique(SI)],           # Line colors
         lwd = 0)                 # Line width
  
  
  
  library(ade4)
  unique(DepthMD[,2])
  ROWNamova<-rownames(DepthMD[DepthMD[,2]=='BCTV-Svr',])
  Structure<-as.data.frame(as.vector(DepthMD[DepthMD[,2]=='BCTV-Svr',4]))
  colnames(Structure)<-'State'
  SAMP<-t(as.data.frame(Depth2[ROWNamova,]))
  colnames(SAMP)<-seq(1,length(SAMP[1,]))
  SAMP1<-cbind(SAMP,AFV[,4])
  SAMP2<-SAMP1[SAMP1[,20]>0,]
  SAMP3<-SAMP2[,-25]
  SAMP4<-as.matrix(SAMP3)
  SAMP5<-c()
  for (i in seq(1,length(SAMP4[,1]))){
    if (length(unique(SAMP4[i,]))>1){
      SAMP5<-rbind(SAMP5,SAMP4[i,])
    }
  }
  
  SMAP5<-data.frame(SAMP5)
  DIST<-dist(t(SAMP5),method='euclidean')
  heatmap(as.matix(SAMP5),Rowv = NA)
  library("pheatmap")
  pheatmap(SAMP5,cluster_rows=F,clustering_method = 'complete',breaks = 10)
  
  
  for (i in seq(1,length(SAMP5[,1]))){
    AMOV<-amova(samples = SAMP5[i,],structures = Structure)
    print(AMOV$componentsofcovariance)
  }
  
  AMOV<-amova(samples=data.frame(SAMP5),structures = Structure)
  AMOV$results
  AMOV$componentsofcovariance
  AMOV$statphi
  
  library(dendextend)
  hc<-hclust(DIST)
  dhc <- as.dendrogram(hc)
  plot(dhc)
  Structure
  dend2 <- color_branches(dhc, k = 6)
  dend2<-set(dend2, "branches_lwd", 3)
  plot(dend2)
  
  # creates a single item vector of the clusters    
  myclusters <- cutree(dhc, k=12)
  
  # make the dataframe of two columns cluster number and label
  clusterDF <-  data.frame(Cluster = as.numeric(unlist(myclusters)),
                           Branch = names(myclusters))
  clusterDF[order(clusterDF[,1]),]
  ROWNamova[clusterDF[order(clusterDF[,1]),2]]
  cbind(DepthMD[ROWNamova[clusterDF[order(clusterDF[,1]),2]],],clusterDF[order(clusterDF[,1]),])
  
  
  # sort by cluster ascending
  clusterDFSort <- clusterDF %>% arrange(Cluster)
  Structure
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  #WORLAND
  {
    unique(DepthMD[,2])
    ROWNamova<-rownames(DepthMD[DepthMD[,2]=='BCTV-Wor',])
    Structure<-as.data.frame(as.vector(DepthMD[DepthMD[,2]=='BCTV-Wor',4]))
    colnames(Structure)<-'State'
    SAMP<-t(as.data.frame(Depth2[ROWNamova,]))
    colnames(SAMP)<-seq(1,length(SAMP[1,]))
    SAMP1<-cbind(SAMP,AFV[,4])
    SAMP2<-SAMP1[SAMP1[,20]>0,]
    SAMP3<-SAMP2[,-25]
    SAMP4<-as.matrix(SAMP)
    SAMP5<-c()
    for (i in seq(1,length(SAMP4[,1]))){
      if (length(unique(SAMP4[i,]))>1){
        SAMP5<-rbind(SAMP5,SAMP4[i,])
      }
    }
    
    SMAP5<-data.frame(SAMP5)
    DIST<-dist(t(SAMP5),method='euclidean')
    #install.packages('dendextend')
    library(dendextend)
    hc<-hclust(DIST)
    dhc <- as.dendrogram(hc)
    dend2 <- color_branches(dhc, k = 5)
    dend2<-set(dend2, "branches_lwd", 3)
    plot(dend2)
    get_leaves_branches_col(dend2)
    # creates a single item vector of the clusters    
    myclusters <- cutree(dhc, k=5)
    # make the dataframe of two columns cluster number and label
    clusterDF <-  data.frame(Cluster = as.numeric(unlist(myclusters)),
                             Branch = names(myclusters))
    clusterDF[order(clusterDF[,1]),]
    ROWNamova[clusterDF[order(clusterDF[,1]),2]]
    CLusterz<-cbind(DepthMD[ROWNamova[clusterDF[order(clusterDF[,1]),2]],],clusterDF[order(clusterDF[,1]),])
    unique(CLusterz[,9])
    
    
    FstALL<-data.frame()
    Cor1<-data.frame()
    AFV<-data.frame(1,2,3,4)
    colnames(AFV)<-c('1','2','3','4')
    clust1<-CLusterz[CLusterz[,9]=='1',]
    clust2<-CLusterz[CLusterz[,9]=='2',]
    clust3<-CLusterz[CLusterz[,9]=='3',]
    clust4<-CLusterz[CLusterz[,9]=='4',]
    clust5<-CLusterz[CLusterz[,9]=='5',]
    clust6<-CLusterz[CLusterz[,9]=='6',]
    
    for (i in seq(1,3175)){
      i<-i
      print(i)
      X<-as.vector(Depth2[rownames(clust1),i])
      Xa<-(length(X[X==1])/length(X))^2
      Xt<-(length(X[X==2])/length(X))^2
      Xg<-(length(X[X==3])/length(X))^2
      Xc<-(length(X[X==4])/length(X))^2
      Xx<-(length(X[X==5])/length(X))^2
      SX<-1-sum(Xa,Xt,Xg,Xc,Xx)
      V1<-c(Xa,Xt,Xg,Xc,Xx)
      SUMX<-(SX)
      
      X1<-as.vector(Depth2[rownames(clust2),i])
      Xa1<-(length(X1[X1==1])/length(X1))^2
      Xt1<-(length(X1[X1==2])/length(X1))^2
      Xg1<-(length(X1[X1==3])/length(X1))^2
      Xc1<-(length(X1[X1==4])/length(X1))^2
      Xx1<-(length(X1[X1==5])/length(X1))^2
      SX1<-1-sum(Xa1,Xt1,Xg1,Xc1,Xx1)
      SUMX1<-(SX1)
      FstSite<-((SX+SX1)-(SX1))/(SX+SX1)
      V2<-c(Xa1,Xt1,Xg1,Xc1,Xx1)
      C1<-cor(V1,V2)
      
      X2<-as.vector(Depth2[rownames(clust3),i])
      Xa2<-(length(X2[X2==1])/length(X2))^2
      Xt2<-(length(X2[X2==2])/length(X2))^2
      Xg2<-(length(X2[X2==3])/length(X2))^2
      Xc2<-(length(X2[X2==4])/length(X2))^2
      Xx2<-(length(X2[X2==5])/length(X2))^2
      SX2<-1-sum(Xa2,Xt2,Xg2,Xc2,Xx2)
      V3<-c(Xa2,Xt2,Xg2,Xc2,Xx2)
      
      X3<-as.vector(Depth2[rownames(clust4),i])
      Xa3<-(length(X3[X3==1])/length(X3))^2
      Xt3<-(length(X3[X3==2])/length(X3))^2
      Xg3<-(length(X3[X3==3])/length(X3))^2
      Xc3<-(length(X3[X3==4])/length(X3))^2
      Xx3<-(length(X3[X3==5])/length(X3))^2
      SX3<-1-sum(Xa3,Xt3,Xg3,Xc3,Xx3)
      V4<-c(Xa3,Xt3,Xg3,Xc3,Xx3)
      
      
      C12id<-cor(V1,V2)
      C13id<-cor(V1,V3)
      C14id<-cor(V1,V4)
      
      
      C21mt<-cor(V2,V1)
      C23mt<-cor(V2,V3)
      C24mt<-cor(V2,V4)
      
      
      C31or<-cor(V3,V1)
      C32or<-cor(V3,V2)
      C34or<-cor(V3,V4)
      
      
      C41wy<-cor(V4,V1)
      C42wy<-cor(V4,V2)
      C43wy<-cor(V4,V3)
      
      
      C1cor<-c(C12id,C13id,C14id,C21mt,C23mt,C24mt,
               C31or,C32or,C34or,C41wy,C42wy,C43wy)
      
      AFV1<-c(SX,SX1,SX2,SX3)
      AFV<-rbind(AFV,AFV1)
      Cor1<-rbind(Cor1,C1cor)
    }
    AFV<-AFV[-1,]
    colMeans(AFV)
    
    Cor1<-(1-abs(Cor1))
    
    
    Cor1
    clust1cor<-rowMeans(Cor1[,c(1,2,3)])
    clust2cor<-rowMeans(Cor1[,c(4,5,6)])
    clust3cor<-rowMeans(Cor1[,c(7,8,9)])
    clust4cor<-rowMeans(Cor1[,c(10,11,12)])
    
    clust1cor[clust1cor<0.6]=0
    clust2cor[clust2cor<0.6]=0
    clust3cor[clust3cor<0.6]=0
    clust4cor[clust4cor<0.6]=0
    
    colMeans(cbind(clust1cor,clust2cor,clust3cor,clust4cor))
    
    
    par(mfrow=c(1,1))
    plot(NULL,ylim=c(0,40),xlim = c(0,3175))
    
    ORFF<-ORFs[ORFs[,1]=='giWOR',c(13,14)]
    
    lorf<-c()
    for (i in seq(1,length(ORFF[,1]))){
      print(max(c(ORFF[i,1],ORFF[i,2]))-min(c(ORFF[i,1],ORFF[i,2])))
      lorf<-c(lorf,max(c(ORFF[i,1],ORFF[i,2]))-min(c(ORFF[i,1],ORFF[i,2])))
    }
    olorf<-order(lorf)
    oolorf<-olorf[seq(1,length(lorf)-7)]
    ORFF<-ORFF[oolorf,]
    
    
    lorf<-c()
    for (p in seq(1,length(ORFF[,1]))){
      ll<-length(seq(ORFF[p,1],ORFF[p,2]))
      lorf<-c(lorf,ll)
    }
    
    lorf<-c()
    for (p in seq(1,length(ORFF[,1]))){
      ll<-(length(seq(ORFF[p,1],ORFF[p,2])))
      lorf<-c(lorf,ll)
      
    }
    lorf
    lorfo<-order(lorf,decreasing = F)
    
    Pdf<-c(rep(0,3175))
    Pfd1<-c(0,0,0)
    for (p in lorfo){
      print(ORFF[p,])
      P1<-c(rep(0,3175))
      P1[seq(ORFF[p,1],ORFF[p,2])]=1
      Pdf<-rbind(Pdf,P1)
      print((colSums(Pdf[,seq(ORFF[p,1],ORFF[p,2])])))
      MCS<-max(colSums(Pdf[,seq(ORFF[p,1],ORFF[p,2])]))
      print(MCS)
      ORFFF<-(cbind(ORFF[p,],MCS))
      Pfd1<-rbind(Pfd1,ORFFF)
    }
    
    Pfd1<-Pfd1[-1,]
    inens<-c()
    for (p in seq(1,length(Pfd1[,1]))){
      inens<-c(inens,1)
    }
    
    
    for (p in seq(1,length(Pfd1[,1]))){
      color_transparent <- adjustcolor('grey', alpha.f = inens[p])
      xbottom <- Pfd1[p,1]
      ybottom <- Pfd1[p,3]+19
      xtop <- Pfd1[p,2]
      ytop<- Pfd1[p,3]+20
      xpol <- c(xbottom,xbottom,xtop,xtop)
      ypol <- c(ybottom,ytop,ytop,ybottom)
      polygon(x=cuts(t(xpol)), y=cuts(t(ypol)), col = color_transparent, border ='black')
    }
    
    setwd('z:/2020/Raj_Paper/CleanData/blast_orfs/')
    #ADD sncRNA
    sncRNA<-read.table('out.snc.blast')
    Pfd1<-sncRNA[,c(9,10)]
    for (p in seq(1,length(Pfd1[,1]))){
      color_transparent <- adjustcolor('black', alpha.f = .5)
      xbottom <- Pfd1[p,1]
      ybottom <- 18
      xtop <- Pfd1[p,2]
      ytop<- 19
      xpol <- c(xbottom,xbottom,xtop,xtop)
      ypol <- c(ybottom,ytop,ytop,ybottom)
      polygon(x=cuts(t(xpol)), y=cuts(t(ypol)), col = color_transparent, border =NA)
    }
    
    
    
    for (p in seq(1,length(clust1cor))){
      inens<- clust1cor[p]*2
      color_transparent <- adjustcolor('orange', alpha.f = inens)
      Gm<-p
      GM<-p+1
      xbottom <- Gm
      ybottom <- 17
      xtop <- GM
      ytop<- 17.5
      xpol <- c(xbottom,xbottom,xtop,xtop)
      ypol <- c(ybottom,ytop,ytop,ybottom)
      polygon(x=cuts(t(xpol)), y=cuts(t(ypol)), col = color_transparent, border =NA)
    }
    
    for (p in seq(1,length(clust2cor))){
      inens<- clust2cor[p]*2
      color_transparent <- adjustcolor('blue', alpha.f = inens)
      Gm<-p
      GM<-p+1
      xbottom <- Gm
      ybottom <- 16.5
      xtop <- GM
      ytop<- 17
      xpol <- c(xbottom,xbottom,xtop,xtop)
      ypol <- c(ybottom,ytop,ytop,ybottom)
      polygon(x=cuts(t(xpol)), y=cuts(t(ypol)), col = color_transparent, border =NA)
    }
    
    for (p in seq(1,length(clust3cor))){
      inens<- clust3cor[p]*2
      color_transparent <- adjustcolor('green', alpha.f = inens)
      Gm<-p
      GM<-p+1
      xbottom <- Gm
      ybottom <- 16
      xtop <- GM
      ytop<- 16.5
      xpol <- c(xbottom,xbottom,xtop,xtop)
      ypol <- c(ybottom,ytop,ytop,ybottom)
      polygon(x=cuts(t(xpol)), y=cuts(t(ypol)), col = color_transparent, border =NA)
    }
    
    for (p in seq(1,length(clust4cor))){
      inens<- clust4cor[p]*2
      color_transparent <- adjustcolor('red', alpha.f = inens)
      Gm<-p
      GM<-p+1
      xbottom <- Gm
      ybottom <- 16
      xtop <- GM
      ytop<- 15.5
      xpol <- c(xbottom,xbottom,xtop,xtop)
      ypol <- c(ybottom,ytop,ytop,ybottom)
      polygon(x=cuts(t(xpol)), y=cuts(t(ypol)), col = color_transparent, border =NA)
    }
    
    #BCTV genes
    {
      color_transparent <- adjustcolor('grey0', alpha.f = 0.5)
      sncRNA<-read.table('BCTVgenes_prot.blast.out')
      Gm<-c(min(sncRNA[sncRNA[,1]=='V1',c(9,10)]))
      GM<-c(max(sncRNA[sncRNA[,1]=='V1',c(9,10)]))
      xbottom <- Gm
      ybottom <- 5
      xtop <- GM
      ytop<- 6
      xpol <- c(xbottom,xbottom,xtop,xtop)
      ypol <- c(ybottom,ytop,ytop,ybottom)
      polygon(x=cuts(t(xpol)), y=cuts(t(ypol)), col = color_transparent, border ='black')
      
      Gm<-c(min(sncRNA[sncRNA[,1]=='V2',c(9,10)]))
      GM<-c(max(sncRNA[sncRNA[,1]=='V2',c(9,10)]))
      xbottom <- Gm
      ybottom <- 6
      xtop <- GM
      ytop<- 7
      xpol <- c(xbottom,xbottom,xtop,xtop)
      ypol <- c(ybottom,ytop,ytop,ybottom)
      polygon(x=cuts(t(xpol)), y=cuts(t(ypol)), col = color_transparent, border ='black')
      
      Gm<-c(min(sncRNA[sncRNA[,1]=='V3',c(9,10)]))
      GM<-c(max(sncRNA[sncRNA[,1]=='V3',c(9,10)]))
      xbottom <- Gm
      ybottom <- 5
      xtop <- GM
      ytop<- 6
      xpol <- c(xbottom,xbottom,xtop,xtop)
      ypol <- c(ybottom,ytop,ytop,ybottom)
      polygon(x=cuts(t(xpol)), y=cuts(t(ypol)), col = color_transparent, border ='black')
      
      color_transparent <- adjustcolor('darkslategray', alpha.f = 0.5)
      Gm<-c(min(sncRNA[sncRNA[,1]=='C1',c(9,10)]))
      GM<-c(max(sncRNA[sncRNA[,1]=='C1',c(9,10)]))
      xbottom <- Gm
      ybottom <- 5
      xtop <- GM
      ytop<- 6
      xpol <- c(xbottom,xbottom,xtop,xtop)
      ypol <- c(ybottom,ytop,ytop,ybottom)
      polygon(x=cuts(t(xpol)), y=cuts(t(ypol)), col = color_transparent, border ='black')
      
      Gm<-c(min(sncRNA[sncRNA[,1]=='C2',c(9,10)]))
      GM<-c(max(sncRNA[sncRNA[,1]=='C2',c(9,10)]))
      xbottom <- Gm
      ybottom <- 6
      xtop <- GM
      ytop<- 7
      xpol <- c(xbottom,xbottom,xtop,xtop)
      ypol <- c(ybottom,ytop,ytop,ybottom)
      polygon(x=cuts(t(xpol)), y=cuts(t(ypol)), col = color_transparent, border ='black')
      
      sncRNA[sncRNA[,1]=='C3',c(9,10)]
      Gm<-c(min(sncRNA[sncRNA[,1]=='C3',c(9,10)]))
      GM<-c(max(sncRNA[sncRNA[,1]=='C3',c(9,10)]))
      xbottom <- Gm
      ybottom <- 5
      xtop <- GM
      ytop<- 6
      xpol <- c(xbottom,xbottom,xtop,xtop)
      ypol <- c(ybottom,ytop,ytop,ybottom)
      polygon(x=cuts(t(xpol)), y=cuts(t(ypol)), col = color_transparent, border ='black')
      
      sncRNA[sncRNA[,1]=='C4',c(9,10)]
      Gm<-c(min(sncRNA[sncRNA[,1]=='C4',c(9,10)]))
      GM<-c(max(sncRNA[sncRNA[,1]=='C4',c(9,10)]))
      xbottom <- Gm
      ybottom <- 6
      xtop <- GM
      ytop<- 7
      xpol <- c(xbottom,xbottom,xtop,xtop)
      ypol <- c(ybottom,ytop,ytop,ybottom)
      polygon(x=cuts(t(xpol)), y=cuts(t(ypol)), col = color_transparent, border ='black')
    }
    
    unique(DepthMD[,2])
    ROWNamova<-CLusterz[,1]
    Structure<-as.data.frame(as.vector(CLusterz[,9]))
    colnames(Structure)<-'Clusters'
    cbind(rep('C',19),Structure)
    stuct<-paste(CLusterz[,9],rep('C',19),sep = '')
    
    ##Reatin Rownames for plots later
    SAMP<-t(as.data.frame(Depth2[ROWNamova,]))
    colnames(SAMP)<-seq(1,length(SAMP[1,]))
    rownames(SAMP)<-seq(1,length(SAMP[,1]))
    SAMP1<-cbind(SAMP,rowMeans(AFV))
    SAMP2<-SAMP1[SAMP1[,20]>0,]
    SAMP3<-SAMP2[,-20]
    SAMP4<-as.matrix(SAMP)
    SAMP5<-c()
    for (i in seq(1,length(SAMP4[,1]))){
      if (length(unique(SAMP4[i,]))>1){
        SAMP5<-rbind(SAMP5,SAMP4[i,])
      }
    }
    
    NCC1t<-data.frame()
    for (i in seq(1,length(SAMP5[,1]))){
      NCC1<-cbind(rep(i,5),c(1,2,3,4,5))
      for (l in seq(1,length(SAMP5[1,]))){
        n1=0
        n2=0
        n3=0
        n4=0
        n5=0
        Nc<-c()
        SAMP7<-SAMP5[i,l]
        n0<-i
        n1<-length(SAMP7[SAMP7=='1'])
        n2<-length(SAMP7[SAMP7=='2'])
        n3<-length(SAMP7[SAMP7=='3'])
        n4<-length(SAMP7[SAMP7=='4'])
        n5<-length(SAMP7[SAMP7=='5'])
        Nc<-c(n1,n2,n3,n4,n5)
        NCC1<-cbind(NCC1,Nc)
      }
      NCC1t<-rbind(NCC1t,NCC1)
      #colnames(NCC1t)<-c('1','2','3','4','5')
    }
    
    rownames(NCC1t)<-paste(NCC1t[,1],NCC1t[,2],sep = '_')
    gsub(pattern = '_',x = rownames(NCC1t),replacement = ' ')
    NCC1t<-NCC1t[,c(-1,-2)]
    
    
    library(ade4)
    AMOVA_out<-data.frame()
    for (i in seq(1,length(NCC1t[,1]))){
      Amovdf<-NCC1t
      AMOV<-amova(samples=data.frame(Amovdf),structures = data.frame(stuct))
      print(AMOV$componentsofcovariance[1,])
      AMOVA_out<-rbind(AMOVA_out,AMOV$componentsofcovariance[1,])
    }
    
    NCC1t<-NCC1t[rowMeans(NCC1t)>0,]
    AMOV<-amova(samples=data.frame(NCC1t),structures = data.frame(stuct))
    
    ###GLM approch
    stuct<-data.frame(stuct)
    stuct<-cbind(rownames(stuct),stuct)
    colsClust<-c(1,2,3,4,5)
    NCC1tP1<-NCC1t[,colsClust]
    NCC1tP2<-NCC1t[,-colsClust]
    
    
    
    TT<-c(1,2,3,4)
    for (i in seq(1,length(NCC1tP2[,1]))){
      Tt<-t.test(NCC1tP1[i,],NCC1tP2[i,])
      Ta<-c(i,Tt$p.value,Tt$estimate)
      TT<-rbind(TT,Ta)
    }
    TTo<-TT[order(TT[,2]),]
    TTo<-TTo[TTo[,2]<0.05,]
    Ortest<-as.vector(TT[order(TT[,2]),][,1])
    NCC1t[TTo[,1],]
    TTo1<-cbind(NCC1t[TTo[,1],],TTo)
    TTo2<-TTo1[TTo1[,22]>0,]
    TTo2<-na.omit(TTo2)
    TTo1[order(TTo1$V1,decreasing = F),]
    vec_pos<-c()
    vec_id<-c()
    for (i in rownames(TTo2)){
      vec_pos<-c(vec_pos,strsplit(i,split = "_")[[1]][1])
      vec_id<-c(vec_id,strsplit(i,split = "_")[[1]][2])
    }
    
    vec_pos<-as.numeric(vec_pos)
    TTo2<-cbind(TTo2,vec_pos,vec_id)
    
    for (p in seq(1,length(TTo2[,1]))){
      color_transparent <- adjustcolor('#917600', alpha.f = 0.7)
      Gm<-TTo2[p,24]
      GM<-TTo2[p,24]+1
      xbottom <- Gm
      ybottom <- 16.5
      xtop <- GM
      ytop<- 17.5
      xpol <- c(xbottom,xbottom,xtop,xtop)
      ypol <- c(ybottom,ytop,ytop,ybottom)
      polygon(x=cuts(t(xpol)), y=cuts(t(ypol)), col = color_transparent, border =NA)
    }
    TTfinalTable1<-TTo2
    #cluster 2
    stuct<-data.frame(stuct)
    stuct<-cbind(rownames(stuct),stuct)
    colsClust<-c(6,7,8,9,10,11)
    NCC1tP1<-NCC1t[,colsClust]
    NCC1tP2<-NCC1t[,-colsClust]
    
    TT<-c(1,2,3,4)
    for (i in seq(1,length(NCC1tP2[,1]))){
      Tt<-t.test(NCC1tP1[i,],NCC1tP2[i,])
      Ta<-c(i,Tt$p.value,Tt$estimate)
      TT<-rbind(TT,Ta)
    }
    TTo<-TT[order(TT[,2]),]
    TTo<-TTo[TTo[,2]<0.05,]
    Ortest<-as.vector(TT[order(TT[,2]),][,1])
    NCC1t[TTo[,1],]
    TTo1<-cbind(NCC1t[TTo[,1],],TTo)
    TTo2<-TTo1[TTo1[,22]>0,]
    TTo2<-na.omit(TTo2)
    TTo1[order(TTo1$V1,decreasing = F),]
    vec_pos<-c()
    vec_id<-c()
    for (i in rownames(TTo2)){
      vec_pos<-c(vec_pos,strsplit(i,split = "_")[[1]][1])
      vec_id<-c(vec_id,strsplit(i,split = "_")[[1]][2])
    }
    
    vec_pos<-as.numeric(vec_pos)
    TTo2<-cbind(TTo2,vec_pos,vec_id)
    
    for (p in seq(1,length(TTo2[,1]))){
      color_transparent <- adjustcolor('#008FB7', alpha.f = 0.7)
      Gm<-TTo2[p,24]
      GM<-TTo2[p,24]+1
      xbottom <- Gm
      ybottom <- 15.5
      xtop <- GM
      ytop<- 16.5
      xpol <- c(xbottom,xbottom,xtop,xtop)
      ypol <- c(ybottom,ytop,ytop,ybottom)
      polygon(x=cuts(t(xpol)), y=cuts(t(ypol)), col = color_transparent, border =NA)
    }
    TTfinalTable2<-TTo2
    
    #cluster 3
    stuct<-data.frame(stuct)
    stuct<-cbind(rownames(stuct),stuct)
    colsClust<-c(11,12)
    NCC1tP1<-NCC1t[,colsClust]
    NCC1tP2<-NCC1t[,-colsClust]
    
    TT<-c(1,2,3,4)
    for (i in seq(1,length(NCC1tP2[,1]))){
      Tt<-t.test(NCC1tP1[i,],NCC1tP2[i,])
      Ta<-c(i,Tt$p.value,Tt$estimate)
      TT<-rbind(TT,Ta)
    }
    TTo<-TT[order(TT[,2]),]
    TTo<-TTo[TTo[,2]<0.05,]
    Ortest<-as.vector(TT[order(TT[,2]),][,1])
    NCC1t[TTo[,1],]
    TTo1<-cbind(NCC1t[TTo[,1],],TTo)
    TTo2<-TTo1[TTo1[,22]>0,]
    TTo2<-na.omit(TTo2)
    TTo1[order(TTo1$V1,decreasing = F),]
    vec_pos<-c()
    vec_id<-c()
    for (i in rownames(TTo2)){
      vec_pos<-c(vec_pos,strsplit(i,split = "_")[[1]][1])
      vec_id<-c(vec_id,strsplit(i,split = "_")[[1]][2])
    }
    
    vec_pos<-as.numeric(vec_pos)
    TTo2<-cbind(TTo2,vec_pos,vec_id)
    
    for (p in seq(1,length(TTo2[,1]))){
      color_transparent <- adjustcolor('#009232', alpha.f = 0.7)
      Gm<-TTo2[p,24]
      GM<-TTo2[p,24]+1
      xbottom <- Gm
      ybottom <- 14.5
      xtop <- GM
      ytop<- 15.5 
      xpol <- c(xbottom,xbottom,xtop,xtop)
      ypol <- c(ybottom,ytop,ytop,ybottom)
      polygon(x=cuts(t(xpol)), y=cuts(t(ypol)), col = color_transparent, border =NA)
    }
    TTfinalTable3<-TTo2
    #cluster 4
    stuct<-data.frame(stuct)
    stuct<-cbind(rownames(stuct),stuct)
    colsClust<-c(13,15,16,17,18)
    NCC1tP1<-NCC1t[,colsClust]
    NCC1tP2<-NCC1t[,-colsClust]
    
    TT<-c(1,2,3,4)
    for (i in seq(1,length(NCC1tP2[,1]))){
      Tt<-t.test(NCC1tP1[i,],NCC1tP2[i,])
      Ta<-c(i,Tt$p.value,Tt$estimate)
      TT<-rbind(TT,Ta)
    }
    TTo<-TT[order(TT[,2]),]
    TTo<-TTo[TTo[,2]<0.05,]
    Ortest<-as.vector(TT[order(TT[,2]),][,1])
    NCC1t[TTo[,1],]
    TTo1<-cbind(NCC1t[TTo[,1],],TTo)
    TTo2<-TTo1[TTo1[,22]>0,]
    TTo2<-na.omit(TTo2)
    TTo1[order(TTo1$V1,decreasing = F),]
    vec_pos<-c()
    vec_id<-c()
    for (i in rownames(TTo2)){
      vec_pos<-c(vec_pos,strsplit(i,split = "_")[[1]][1])
      vec_id<-c(vec_id,strsplit(i,split = "_")[[1]][2])
    }
    
    vec_pos<-as.numeric(vec_pos)
    TTo2<-cbind(TTo2,vec_pos,vec_id)
    
    for (p in seq(1,length(TTo2[,1]))){
      color_transparent <- adjustcolor('#A352D1', alpha.f = 0.7)
      Gm<-TTo2[p,24]
      GM<-TTo2[p,24]+1
      xbottom <- Gm
      ybottom <- 14.5
      xtop <- GM
      ytop<- 13.5
      xpol <- c(xbottom,xbottom,xtop,xtop)
      ypol <- c(ybottom,ytop,ytop,ybottom)
      polygon(x=cuts(t(xpol)), y=cuts(t(ypol)), col = color_transparent, border =NA)
    }
    TTfinalTable4<-TTo2
    
  }    
    

  #CALOGAN
  {
    unique(DepthMD[,2])
    ROWNamova<-rownames(DepthMD[DepthMD[,2]=='BCTV-CALogan',])
    Structure<-as.data.frame(as.vector(DepthMD[DepthMD[,2]=='BCTV-CALogan',4]))
    colnames(Structure)<-'State'
    SAMP<-t(as.data.frame(Depth2[ROWNamova,]))
    colnames(SAMP)<-seq(1,length(SAMP[1,]))
    SAMP1<-cbind(SAMP,AFV[,4])
    SAMP2<-SAMP1[SAMP1[,14]>0,]
    SAMP3<-SAMP2[,-14]
    SAMP4<-as.matrix(SAMP)
    SAMP5<-c()
    for (i in seq(1,length(SAMP4[,1]))){
      if (length(unique(SAMP4[i,]))>1){
        SAMP5<-rbind(SAMP5,SAMP4[i,])
      }
    }
    
    SMAP5<-data.frame(SAMP5)
    DIST<-dist(t(SAMP5),method='euclidean')
    #install.packages('dendextend')
    library(dendextend)
    hc<-hclust(DIST)
    dhc <- as.dendrogram(hc)
    dend2 <- color_branches(dhc, k = 4)
    dend2<-set(dend2, "branches_lwd", 3)
    plot(dend2)
    get_leaves_branches_col(dend2)
    # creates a single item vector of the clusters    
    myclusters <- cutree(dhc, k=4)
    # make the dataframe of two columns cluster number and label
    clusterDF <-  data.frame(Cluster = as.numeric(unlist(myclusters)),
                             Branch = names(myclusters))
    clusterDF[order(clusterDF[,1]),]
    ROWNamova[clusterDF[order(clusterDF[,1]),2]]
    CLusterz<-cbind(DepthMD[ROWNamova[clusterDF[order(clusterDF[,1]),2]],],clusterDF[order(clusterDF[,1]),])
    unique(CLusterz[,9])
    
    
    FstALL<-data.frame()
    Cor1<-data.frame()
    AFV<-data.frame(1,2,3,4)
    colnames(AFV)<-c('1','2','3','4')
    clust1<-CLusterz[CLusterz[,9]=='1',]
    clust2<-CLusterz[CLusterz[,9]=='2',]
    clust3<-CLusterz[CLusterz[,9]=='3',]
    clust4<-CLusterz[CLusterz[,9]=='4',]

    
    for (i in seq(1,3175)){
      i<-i
      print(i)
      X<-as.vector(Depth2[rownames(clust1),i])
      Xa<-(length(X[X==1])/length(X))^2
      Xt<-(length(X[X==2])/length(X))^2
      Xg<-(length(X[X==3])/length(X))^2
      Xc<-(length(X[X==4])/length(X))^2
      Xx<-(length(X[X==5])/length(X))^2
      SX<-1-sum(Xa,Xt,Xg,Xc,Xx)
      V1<-c(Xa,Xt,Xg,Xc,Xx)
      SUMX<-(SX)
      
      X1<-as.vector(Depth2[rownames(clust2),i])
      Xa1<-(length(X1[X1==1])/length(X1))^2
      Xt1<-(length(X1[X1==2])/length(X1))^2
      Xg1<-(length(X1[X1==3])/length(X1))^2
      Xc1<-(length(X1[X1==4])/length(X1))^2
      Xx1<-(length(X1[X1==5])/length(X1))^2
      SX1<-1-sum(Xa1,Xt1,Xg1,Xc1,Xx1)
      SUMX1<-(SX1)
      FstSite<-((SX+SX1)-(SX1))/(SX+SX1)
      V2<-c(Xa1,Xt1,Xg1,Xc1,Xx1)
      C1<-cor(V1,V2)
      
      X2<-as.vector(Depth2[rownames(clust3),i])
      Xa2<-(length(X2[X2==1])/length(X2))^2
      Xt2<-(length(X2[X2==2])/length(X2))^2
      Xg2<-(length(X2[X2==3])/length(X2))^2
      Xc2<-(length(X2[X2==4])/length(X2))^2
      Xx2<-(length(X2[X2==5])/length(X2))^2
      SX2<-1-sum(Xa2,Xt2,Xg2,Xc2,Xx2)
      V3<-c(Xa2,Xt2,Xg2,Xc2,Xx2)
      
      X3<-as.vector(Depth2[rownames(clust4),i])
      Xa3<-(length(X3[X3==1])/length(X3))^2
      Xt3<-(length(X3[X3==2])/length(X3))^2
      Xg3<-(length(X3[X3==3])/length(X3))^2
      Xc3<-(length(X3[X3==4])/length(X3))^2
      Xx3<-(length(X3[X3==5])/length(X3))^2
      SX3<-1-sum(Xa3,Xt3,Xg3,Xc3,Xx3)
      V4<-c(Xa3,Xt3,Xg3,Xc3,Xx3)
      
      
      C12id<-cor(V1,V2)
      C13id<-cor(V1,V3)
      C14id<-cor(V1,V4)
      
      
      C21mt<-cor(V2,V1)
      C23mt<-cor(V2,V3)
      C24mt<-cor(V2,V4)
      
      
      C31or<-cor(V3,V1)
      C32or<-cor(V3,V2)
      C34or<-cor(V3,V4)
      
      
      C41wy<-cor(V4,V1)
      C42wy<-cor(V4,V2)
      C43wy<-cor(V4,V3)
      
      
      C1cor<-c(C12id,C13id,C14id,C21mt,C23mt,C24mt,
               C31or,C32or,C34or,C41wy,C42wy,C43wy)
      
      AFV1<-c(SX,SX1,SX2,SX3)
      AFV<-rbind(AFV,AFV1)
      Cor1<-rbind(Cor1,C1cor)
    }
    AFV<-AFV[-1,]
    colMeans(AFV)
    
    Cor1<-(1-abs(Cor1))
    
    
    Cor1
    clust1cor<-rowMeans(Cor1[,c(1,2,3)])
    clust2cor<-rowMeans(Cor1[,c(4,5,6)])
    clust3cor<-rowMeans(Cor1[,c(7,8,9)])
    clust4cor<-rowMeans(Cor1[,c(10,11,12)])
    
    colMeans(cbind(clust1cor,clust2cor,clust3cor,clust4cor))
    
    clust1cor[clust1cor<0.6]=0
    clust2cor[clust2cor<0.6]=0
    clust3cor[clust3cor<0.6]=0
    clust4cor[clust4cor<0.6]=0
    
    par(mfrow=c(1,1))
    plot(NULL,ylim=c(0,40),xlim = c(0,3175))
    
    ORFF<-ORFs[ORFs[,1]=='giCAL',c(13,14)]
    
    lorf<-c()
    for (i in seq(1,length(ORFF[,1]))){
      print(max(c(ORFF[i,1],ORFF[i,2]))-min(c(ORFF[i,1],ORFF[i,2])))
      lorf<-c(lorf,max(c(ORFF[i,1],ORFF[i,2]))-min(c(ORFF[i,1],ORFF[i,2])))
    }
    olorf<-order(lorf)
    oolorf<-olorf[seq(1,length(lorf)-7)]
    ORFF<-ORFF[oolorf,]
    
    
    lorf<-c()
    for (p in seq(1,length(ORFF[,1]))){
      ll<-length(seq(ORFF[p,1],ORFF[p,2]))
      lorf<-c(lorf,ll)
    }
    
    lorf<-c()
    for (p in seq(1,length(ORFF[,1]))){
      ll<-(length(seq(ORFF[p,1],ORFF[p,2])))
      lorf<-c(lorf,ll)
      
    }
    lorf
    lorfo<-order(lorf,decreasing = F)
    
    Pdf<-c(rep(0,3175))
    Pfd1<-c(0,0,0)
    for (p in lorfo){
      print(ORFF[p,])
      P1<-c(rep(0,3175))
      P1[seq(ORFF[p,1],ORFF[p,2])]=1
      Pdf<-rbind(Pdf,P1)
      print((colSums(Pdf[,seq(ORFF[p,1],ORFF[p,2])])))
      MCS<-max(colSums(Pdf[,seq(ORFF[p,1],ORFF[p,2])]))
      print(MCS)
      ORFFF<-(cbind(ORFF[p,],MCS))
      Pfd1<-rbind(Pfd1,ORFFF)
    }
    
    Pfd1<-Pfd1[-1,]
    inens<-c()
    for (p in seq(1,length(Pfd1[,1]))){
      inens<-c(inens,1)
    }
    
    
    for (p in seq(1,length(Pfd1[,1]))){
      color_transparent <- adjustcolor('grey', alpha.f = inens[p])
      xbottom <- Pfd1[p,1]
      ybottom <- Pfd1[p,3]+19
      xtop <- Pfd1[p,2]
      ytop<- Pfd1[p,3]+20
      xpol <- c(xbottom,xbottom,xtop,xtop)
      ypol <- c(ybottom,ytop,ytop,ybottom)
      polygon(x=cuts(t(xpol)), y=cuts(t(ypol)), col = color_transparent, border ='black')
    }
    
    
    #ADD sncRNA
    sncRNA<-read.table('out.snc.blast')
    Pfd1<-sncRNA[,c(9,10)]
    for (p in seq(1,length(Pfd1[,1]))){
      color_transparent <- adjustcolor('black', alpha.f = .5)
      xbottom <- Pfd1[p,1]
      ybottom <- 18
      xtop <- Pfd1[p,2]
      ytop<- 19
      xpol <- c(xbottom,xbottom,xtop,xtop)
      ypol <- c(ybottom,ytop,ytop,ybottom)
      polygon(x=cuts(t(xpol)), y=cuts(t(ypol)), col = color_transparent, border =NA)
    }
    
#FST    
{   
    for (p in seq(1,length(clust1cor))){
      inens<- clust1cor[p]*1
      color_transparent <- adjustcolor('orange', alpha.f = inens)
      Gm<-p
      GM<-p+1
      xbottom <- Gm
      ybottom <- 17
      xtop <- GM
      ytop<- 17.5
      xpol <- c(xbottom,xbottom,xtop,xtop)
      ypol <- c(ybottom,ytop,ytop,ybottom)
      polygon(x=cuts(t(xpol)), y=cuts(t(ypol)), col = color_transparent, border =NA)
    }
    
    for (p in seq(1,length(clust2cor))){
      inens<- clust2cor[p]*1
      color_transparent <- adjustcolor('blue', alpha.f = inens)
      Gm<-p
      GM<-p+1
      xbottom <- Gm
      ybottom <- 16.5
      xtop <- GM
      ytop<- 17
      xpol <- c(xbottom,xbottom,xtop,xtop)
      ypol <- c(ybottom,ytop,ytop,ybottom)
      polygon(x=cuts(t(xpol)), y=cuts(t(ypol)), col = color_transparent, border =NA)
    }
    
    for (p in seq(1,length(clust3cor))){
      inens<- clust3cor[p]*1
      color_transparent <- adjustcolor('green', alpha.f = inens)
      Gm<-p
      GM<-p+1
      xbottom <- Gm
      ybottom <- 16
      xtop <- GM
      ytop<- 16.5
      xpol <- c(xbottom,xbottom,xtop,xtop)
      ypol <- c(ybottom,ytop,ytop,ybottom)
      polygon(x=cuts(t(xpol)), y=cuts(t(ypol)), col = color_transparent, border =NA)
    }
    
    for (p in seq(1,length(clust4cor))){
      inens<- clust4cor[p]*1
      color_transparent <- adjustcolor('red', alpha.f = inens)
      Gm<-p
      GM<-p+1
      xbottom <- Gm
      ybottom <- 16
      xtop <- GM
      ytop<- 15.5
      xpol <- c(xbottom,xbottom,xtop,xtop)
      ypol <- c(ybottom,ytop,ytop,ybottom)
      polygon(x=cuts(t(xpol)), y=cuts(t(ypol)), col = color_transparent, border =NA)
    }
    
}   
    #BCTV genes
    {
      color_transparent <- adjustcolor('grey0', alpha.f = 0.5)
      sncRNA<-read.table('BCTVgenes_prot.blast.out')
      Gm<-c(min(sncRNA[sncRNA[,1]=='V1',c(9,10)]))
      GM<-c(max(sncRNA[sncRNA[,1]=='V1',c(9,10)]))
      xbottom <- Gm
      ybottom <- 5
      xtop <- GM
      ytop<- 6
      xpol <- c(xbottom,xbottom,xtop,xtop)
      ypol <- c(ybottom,ytop,ytop,ybottom)
      polygon(x=cuts(t(xpol)), y=cuts(t(ypol)), col = color_transparent, border ='black')
      
      Gm<-c(min(sncRNA[sncRNA[,1]=='V2',c(9,10)]))
      GM<-c(max(sncRNA[sncRNA[,1]=='V2',c(9,10)]))
      xbottom <- Gm
      ybottom <- 6
      xtop <- GM
      ytop<- 7
      xpol <- c(xbottom,xbottom,xtop,xtop)
      ypol <- c(ybottom,ytop,ytop,ybottom)
      polygon(x=cuts(t(xpol)), y=cuts(t(ypol)), col = color_transparent, border ='black')
      
      Gm<-c(min(sncRNA[sncRNA[,1]=='V3',c(9,10)]))
      GM<-c(max(sncRNA[sncRNA[,1]=='V3',c(9,10)]))
      xbottom <- Gm
      ybottom <- 5
      xtop <- GM
      ytop<- 6
      xpol <- c(xbottom,xbottom,xtop,xtop)
      ypol <- c(ybottom,ytop,ytop,ybottom)
      polygon(x=cuts(t(xpol)), y=cuts(t(ypol)), col = color_transparent, border ='black')
      
      color_transparent <- adjustcolor('darkslategray', alpha.f = 0.5)
      Gm<-c(min(sncRNA[sncRNA[,1]=='C1',c(9,10)]))
      GM<-c(max(sncRNA[sncRNA[,1]=='C1',c(9,10)]))
      xbottom <- Gm
      ybottom <- 5
      xtop <- GM
      ytop<- 6
      xpol <- c(xbottom,xbottom,xtop,xtop)
      ypol <- c(ybottom,ytop,ytop,ybottom)
      polygon(x=cuts(t(xpol)), y=cuts(t(ypol)), col = color_transparent, border ='black')
      
      Gm<-c(min(sncRNA[sncRNA[,1]=='C2',c(9,10)]))
      GM<-c(max(sncRNA[sncRNA[,1]=='C2',c(9,10)]))
      xbottom <- Gm
      ybottom <- 6
      xtop <- GM
      ytop<- 7
      xpol <- c(xbottom,xbottom,xtop,xtop)
      ypol <- c(ybottom,ytop,ytop,ybottom)
      polygon(x=cuts(t(xpol)), y=cuts(t(ypol)), col = color_transparent, border ='black')
      
      sncRNA[sncRNA[,1]=='C3',c(9,10)]
      Gm<-c(min(sncRNA[sncRNA[,1]=='C3',c(9,10)]))
      GM<-c(max(sncRNA[sncRNA[,1]=='C3',c(9,10)]))
      xbottom <- Gm
      ybottom <- 5
      xtop <- GM
      ytop<- 6
      xpol <- c(xbottom,xbottom,xtop,xtop)
      ypol <- c(ybottom,ytop,ytop,ybottom)
      polygon(x=cuts(t(xpol)), y=cuts(t(ypol)), col = color_transparent, border ='black')
      
      sncRNA[sncRNA[,1]=='C4',c(9,10)]
      Gm<-c(min(sncRNA[sncRNA[,1]=='C4',c(9,10)]))
      GM<-c(max(sncRNA[sncRNA[,1]=='C4',c(9,10)]))
      xbottom <- Gm
      ybottom <- 6
      xtop <- GM
      ytop<- 7
      xpol <- c(xbottom,xbottom,xtop,xtop)
      ypol <- c(ybottom,ytop,ytop,ybottom)
      polygon(x=cuts(t(xpol)), y=cuts(t(ypol)), col = color_transparent, border ='black')
    }
    
    unique(DepthMD[,2])
    ROWNamova<-CLusterz[,1]
    Structure<-as.data.frame(as.vector(CLusterz[,9]))
    colnames(Structure)<-'Clusters'
    cbind(rep('C',19),Structure)
    stuct<-paste(CLusterz[,9],rep('C',19),sep = '')
    
    ##Reatin Rownames for plots later
    SAMP<-t(as.data.frame(Depth2[ROWNamova,]))
    colnames(SAMP)<-seq(1,length(SAMP[1,]))
    rownames(SAMP)<-seq(1,length(SAMP[,1]))
    SAMP1<-cbind(SAMP,rowMeans(AFV))
    SAMP2<-SAMP1[SAMP1[,14]>0,]
    SAMP3<-SAMP2[,-14]
    SAMP4<-as.matrix(SAMP)
    SAMP5<-c()
    for (i in seq(1,length(SAMP4[,1]))){
      if (length(unique(SAMP4[i,]))>1){
        SAMP5<-rbind(SAMP5,SAMP4[i,])
      }
    }
    
    NCC1t<-data.frame()
    for (i in seq(1,length(SAMP5[,1]))){
      NCC1<-cbind(rep(i,5),c(1,2,3,4,5))
      for (l in seq(1,length(SAMP5[1,]))){
        n1=0
        n2=0
        n3=0
        n4=0
        n5=0
        Nc<-c()
        SAMP7<-SAMP5[i,l]
        n0<-i
        n1<-length(SAMP7[SAMP7=='1'])
        n2<-length(SAMP7[SAMP7=='2'])
        n3<-length(SAMP7[SAMP7=='3'])
        n4<-length(SAMP7[SAMP7=='4'])
        n5<-length(SAMP7[SAMP7=='5'])
        Nc<-c(n1,n2,n3,n4,n5)
        NCC1<-cbind(NCC1,Nc)
      }
      NCC1t<-rbind(NCC1t,NCC1)
      #colnames(NCC1t)<-c('1','2','3','4','5')
    }
    
    rownames(NCC1t)<-paste(NCC1t[,1],NCC1t[,2],sep = '_')
    gsub(pattern = '_',x = rownames(NCC1t),replacement = ' ')
    NCC1t<-NCC1t[,c(-1,-2)]
    
    ###GLM approch
    stuct<-data.frame(stuct)
    stuct<-cbind(rownames(stuct),stuct)
    colsClust<-c(1)
    NCC1tP1<-NCC1t[,colsClust]
    NCC1tP2<-NCC1t[,-colsClust]
    
    
    
    TT<-c(1,2,3,4)
    for (i in seq(1,length(NCC1tP2[,1]))){
      Tt<-t.test(NCC1tP1[i,],NCC1tP2[i,])
      Ta<-c(i,Tt$p.value,Tt$estimate)
      TT<-rbind(TT,Ta)
    }
    TTo<-TT[order(TT[,2]),]
    TTo<-TTo[TTo[,2]<0.05,]
    Ortest<-as.vector(TT[order(TT[,2]),][,1])
    NCC1t[TTo[,1],]
    TTo1<-cbind(NCC1t[TTo[,1],],TTo)
    TTo2<-TTo1[TTo1[,16]>0,]
    TTo2<-na.omit(TTo2)
    TTo1[order(TTo1$V1,decreasing = F),]
    vec_pos<-c()
    vec_id<-c()
    for (i in rownames(TTo2)){
      vec_pos<-c(vec_pos,strsplit(i,split = "_")[[1]][1])
      vec_id<-c(vec_id,strsplit(i,split = "_")[[1]][2])
    }
    
    vec_pos<-as.numeric(vec_pos)
    TTo2<-cbind(TTo2,vec_pos,vec_id)
    
    for (p in seq(1,length(TTo2[,1]))){
      color_transparent <- adjustcolor('#CC476B', alpha.f = 1)
      Gm<-TTo2[p,18]
      GM<-TTo2[p,18]+1
      xbottom <- Gm
      ybottom <- 16.5
      xtop <- GM
      ytop<- 17.5
      xpol <- c(xbottom,xbottom,xtop,xtop)
      ypol <- c(ybottom,ytop,ytop,ybottom)
      polygon(x=cuts(t(xpol)), y=cuts(t(ypol)), col = color_transparent, border =NA)
    }
    TTfinalTable1<-TTo2
    #cluster 2
    stuct<-data.frame(stuct)
    stuct<-cbind(rownames(stuct),stuct)
    colsClust<-c(2,3,4,5)
    NCC1tP1<-NCC1t[,colsClust]
    NCC1tP2<-NCC1t[,-colsClust]
    
    TT<-c(1,2,3,4)
    for (i in seq(1,length(NCC1tP2[,1]))){
      Tt<-t.test(NCC1tP1[i,],NCC1tP2[i,])
      Ta<-c(i,Tt$p.value,Tt$estimate)
      TT<-rbind(TT,Ta)
    }
    TTo<-TT[order(TT[,2]),]
    TTo<-TTo[TTo[,2]<0.05,]
    Ortest<-as.vector(TT[order(TT[,2]),][,1])
    NCC1t[TTo[,1],]
    TTo1<-cbind(NCC1t[TTo[,1],],TTo)
    TTo2<-TTo1[TTo1[,16]>0,]
    TTo2<-na.omit(TTo2)
    TTo1[order(TTo1$V1,decreasing = F),]
    vec_pos<-c()
    vec_id<-c()
    for (i in rownames(TTo2)){
      vec_pos<-c(vec_pos,strsplit(i,split = "_")[[1]][1])
      vec_id<-c(vec_id,strsplit(i,split = "_")[[1]][2])
    }
    
    vec_pos<-as.numeric(vec_pos)
    TTo2<-cbind(TTo2,vec_pos,vec_id)
    
    for (p in seq(1,length(TTo2[,1]))){
      color_transparent <- adjustcolor('#009681', alpha.f = 0.7)
      Gm<-TTo2[p,18]
      GM<-TTo2[p,18]+1
      xbottom <- Gm
      ybottom <- 15.5
      xtop <- GM
      ytop<- 16.5
      xpol <- c(xbottom,xbottom,xtop,xtop)
      ypol <- c(ybottom,ytop,ytop,ybottom)
      polygon(x=cuts(t(xpol)), y=cuts(t(ypol)), col = color_transparent, border =NA)
    }
    TTfinalTable2<-TTo2
    
    #cluster 3
    stuct<-data.frame(stuct)
    stuct<-cbind(rownames(stuct),stuct)
    colsClust<-c(6,7)
    NCC1tP1<-NCC1t[,colsClust]
    NCC1tP2<-NCC1t[,-colsClust]
    
    TT<-c(1,2,3,4)
    for (i in seq(1,length(NCC1tP2[,1]))){
      Tt<-t.test(NCC1tP1[i,],NCC1tP2[i,])
      Ta<-c(i,Tt$p.value,Tt$estimate)
      TT<-rbind(TT,Ta)
    }
    TTo<-TT[order(TT[,2]),]
    TTo<-TTo[TTo[,2]<0.05,]
    Ortest<-as.vector(TT[order(TT[,2]),][,1])
    NCC1t[TTo[,1],]
    TTo1<-cbind(NCC1t[TTo[,1],],TTo)
    TTo2<-TTo1[TTo1[,16]>0,]
    TTo2<-na.omit(TTo2)
    TTo1[order(TTo1$V1,decreasing = F),]
    vec_pos<-c()
    vec_id<-c()
    for (i in rownames(TTo2)){
      vec_pos<-c(vec_pos,strsplit(i,split = "_")[[1]][1])
      vec_id<-c(vec_id,strsplit(i,split = "_")[[1]][2])
    }
    
    vec_pos<-as.numeric(vec_pos)
    TTo2<-cbind(TTo2,vec_pos,vec_id)
    
    for (p in seq(1,length(TTo2[,1]))){
      color_transparent <- adjustcolor('#7866D8', alpha.f = 0.7)
      Gm<-TTo2[p,18]
      GM<-TTo2[p,18]+1
      xbottom <- Gm
      ybottom <- 14.5
      xtop <- GM
      ytop<- 15.5 
      xpol <- c(xbottom,xbottom,xtop,xtop)
      ypol <- c(ybottom,ytop,ytop,ybottom)
      polygon(x=cuts(t(xpol)), y=cuts(t(ypol)), col = color_transparent, border =NA)
    }
    TTfinalTable3<-TTo2
    #cluster 4
    stuct<-data.frame(stuct)
    stuct<-cbind(rownames(stuct),stuct)
    colsClust<-c(8,9,10,11,12,13)
    NCC1tP1<-NCC1t[,colsClust]
    NCC1tP2<-NCC1t[,-colsClust]
    
    TT<-c(1,2,3,4)
    for (i in seq(1,length(NCC1tP2[,1]))){
      Tt<-t.test(NCC1tP1[i,],NCC1tP2[i,])
      Ta<-c(i,Tt$p.value,Tt$estimate)
      TT<-rbind(TT,Ta)
    }
    TTo<-TT[order(TT[,2]),]
    TTo<-TTo[TTo[,2]<0.05,]
    Ortest<-as.vector(TT[order(TT[,2]),][,1])
    NCC1t[TTo[,1],]
    TTo1<-cbind(NCC1t[TTo[,1],],TTo)
    TTo2<-TTo1[TTo1[,16]>0,]
    TTo2<-na.omit(TTo2)
    TTo1[order(TTo1$V1,decreasing = F),]
    vec_pos<-c()
    vec_id<-c()
    for (i in rownames(TTo2)){
      vec_pos<-c(vec_pos,strsplit(i,split = "_")[[1]][1])
      vec_id<-c(vec_id,strsplit(i,split = "_")[[1]][2])
    }
    
    vec_pos<-as.numeric(vec_pos)
    TTo2<-cbind(TTo2,vec_pos,vec_id)
    
    for (p in seq(1,length(TTo2[,1]))){
      color_transparent <- adjustcolor('#767F00', alpha.f = 0.7)
      Gm<-TTo2[p,18]
      GM<-TTo2[p,18]+1
      xbottom <- Gm
      ybottom <- 13.5
      xtop <- GM
      ytop<- 14.5
      xpol <- c(xbottom,xbottom,xtop,xtop)
      ypol <- c(ybottom,ytop,ytop,ybottom)
      polygon(x=cuts(t(xpol)), y=cuts(t(ypol)), col = color_transparent, border =NA)
    }
    TTfinalTable4<-TTo2
    
  }    
    
  ###Change Colnumbers for >0 and change colnumbers for plots. 
  
  #Severe
  {
    unique(DepthMD[,2])
    ROWNamova<-rownames(DepthMD[DepthMD[,2]=='BCTV-Svr',])
    Structure<-as.data.frame(as.vector(DepthMD[DepthMD[,2]=='BCTV-Svr',4]))
    colnames(Structure)<-'State'
    SAMP<-t(as.data.frame(Depth2[ROWNamova,]))
    colnames(SAMP)<-seq(1,length(SAMP[1,]))
    SAMP1<-cbind(SAMP,AFV[,4])
    SAMP2<-SAMP1[SAMP1[,20]>0,]
    SAMP3<-SAMP2[,-20]
    SAMP4<-as.matrix(SAMP)
    SAMP5<-c()
    for (i in seq(1,length(SAMP4[,1]))){
      if (length(unique(SAMP4[i,]))>1){
        SAMP5<-rbind(SAMP5,SAMP4[i,])
      }
    }
    
    SMAP5<-data.frame(SAMP5)
    DIST<-dist(t(SAMP5),method='euclidean')
    #install.packages('dendextend')
    library(dendextend)
    hc<-hclust(DIST)
    dhc <- as.dendrogram(hc)
    dend2 <- color_branches(dhc, k = 6)
    dend2<-set(dend2, "branches_lwd", 3)
    plot(dend2)
    get_leaves_branches_col(dend2)
    # creates a single item vector of the clusters    
    myclusters <- cutree(dhc, k=6)
    # make the dataframe of two columns cluster number and label
    clusterDF <-  data.frame(Cluster = as.numeric(unlist(myclusters)),
                             Branch = names(myclusters))
    clusterDF[order(clusterDF[,1]),]
    ROWNamova[clusterDF[order(clusterDF[,1]),2]]
    CLusterz<-cbind(DepthMD[ROWNamova[clusterDF[order(clusterDF[,1]),2]],],clusterDF[order(clusterDF[,1]),])
    unique(CLusterz[,9])
    
    
    FstALL<-data.frame()
    Cor1<-data.frame()
    AFV<-data.frame(1,2,3,4)
    colnames(AFV)<-c('1','2','3','4')
    clust1<-CLusterz[CLusterz[,9]=='1',]
    clust2<-CLusterz[CLusterz[,9]=='2',]
    clust3<-CLusterz[CLusterz[,9]=='4',]
    clust4<-CLusterz[CLusterz[,9]=='5',]
    clust5<-CLusterz[CLusterz[,9]=='5',]
    clust6<-CLusterz[CLusterz[,9]=='6',]
    
    for (i in seq(1,3175)){
      i<-i
      print(i)
      X<-as.vector(Depth2[rownames(clust1),i])
      Xa<-(length(X[X==1])/length(X))^2
      Xt<-(length(X[X==2])/length(X))^2
      Xg<-(length(X[X==3])/length(X))^2
      Xc<-(length(X[X==4])/length(X))^2
      Xx<-(length(X[X==5])/length(X))^2
      SX<-1-sum(Xa,Xt,Xg,Xc,Xx)
      V1<-c(Xa,Xt,Xg,Xc,Xx)
      SUMX<-(SX)
      
      X1<-as.vector(Depth2[rownames(clust2),i])
      Xa1<-(length(X1[X1==1])/length(X1))^2
      Xt1<-(length(X1[X1==2])/length(X1))^2
      Xg1<-(length(X1[X1==3])/length(X1))^2
      Xc1<-(length(X1[X1==4])/length(X1))^2
      Xx1<-(length(X1[X1==5])/length(X1))^2
      SX1<-1-sum(Xa1,Xt1,Xg1,Xc1,Xx1)
      SUMX1<-(SX1)
      FstSite<-((SX+SX1)-(SX1))/(SX+SX1)
      V2<-c(Xa1,Xt1,Xg1,Xc1,Xx1)
      C1<-cor(V1,V2)
      
      X2<-as.vector(Depth2[rownames(clust3),i])
      Xa2<-(length(X2[X2==1])/length(X2))^2
      Xt2<-(length(X2[X2==2])/length(X2))^2
      Xg2<-(length(X2[X2==3])/length(X2))^2
      Xc2<-(length(X2[X2==4])/length(X2))^2
      Xx2<-(length(X2[X2==5])/length(X2))^2
      SX2<-1-sum(Xa2,Xt2,Xg2,Xc2,Xx2)
      V3<-c(Xa2,Xt2,Xg2,Xc2,Xx2)
      
      X3<-as.vector(Depth2[rownames(clust4),i])
      Xa3<-(length(X3[X3==1])/length(X3))^2
      Xt3<-(length(X3[X3==2])/length(X3))^2
      Xg3<-(length(X3[X3==3])/length(X3))^2
      Xc3<-(length(X3[X3==4])/length(X3))^2
      Xx3<-(length(X3[X3==5])/length(X3))^2
      SX3<-1-sum(Xa3,Xt3,Xg3,Xc3,Xx3)
      V4<-c(Xa3,Xt3,Xg3,Xc3,Xx3)
      
      
      C12id<-cor(V1,V2)
      C13id<-cor(V1,V3)
      C14id<-cor(V1,V4)
      
      
      C21mt<-cor(V2,V1)
      C23mt<-cor(V2,V3)
      C24mt<-cor(V2,V4)
      
      
      C31or<-cor(V3,V1)
      C32or<-cor(V3,V2)
      C34or<-cor(V3,V4)
      
      
      C41wy<-cor(V4,V1)
      C42wy<-cor(V4,V2)
      C43wy<-cor(V4,V3)
      
      
      C1cor<-c(C12id,C13id,C14id,C21mt,C23mt,C24mt,
               C31or,C32or,C34or,C41wy,C42wy,C43wy)
      
      AFV1<-c(SX,SX1,SX2,SX3)
      AFV<-rbind(AFV,AFV1)
      Cor1<-rbind(Cor1,C1cor)
    }
    AFV<-AFV[-1,]
    colMeans(AFV)
    
    Cor1<-(1-abs(Cor1))
    
    
    Cor1
    clust1cor<-rowMeans(Cor1[,c(1,2,3)])
    clust2cor<-rowMeans(Cor1[,c(4,5,6)])
    clust3cor<-rowMeans(Cor1[,c(7,8,9)])
    clust4cor<-rowMeans(Cor1[,c(10,11,12)])
    
    colMeans(cbind(clust1cor,clust2cor,clust3cor,clust4cor))
    
    clust1cor[clust1cor<0.6]=0
    clust2cor[clust2cor<0.6]=0
    clust3cor[clust3cor<0.6]=0
    clust4cor[clust4cor<0.6]=0
    
    
    
    
    par(mfrow=c(1,1))
    plot(NULL,ylim=c(0,40),xlim = c(0,3175))
    
    ORFF<-ORFs[ORFs[,1]=='giSEV',c(13,14)]
    
    lorf<-c()
    for (i in seq(1,length(ORFF[,1]))){
      print(max(c(ORFF[i,1],ORFF[i,2]))-min(c(ORFF[i,1],ORFF[i,2])))
      lorf<-c(lorf,max(c(ORFF[i,1],ORFF[i,2]))-min(c(ORFF[i,1],ORFF[i,2])))
    }
    olorf<-order(lorf)
    oolorf<-olorf[seq(1,length(lorf)-7)]
    ORFF<-ORFF[oolorf,]
    
    
    lorf<-c()
    for (p in seq(1,length(ORFF[,1]))){
      ll<-length(seq(ORFF[p,1],ORFF[p,2]))
      lorf<-c(lorf,ll)
    }
    
    lorf<-c()
    for (p in seq(1,length(ORFF[,1]))){
      ll<-(length(seq(ORFF[p,1],ORFF[p,2])))
      lorf<-c(lorf,ll)
      
    }
    lorf
    lorfo<-order(lorf,decreasing = F)
    
    Pdf<-c(rep(0,3175))
    Pfd1<-c(0,0,0)
    for (p in lorfo){
      print(ORFF[p,])
      P1<-c(rep(0,3175))
      P1[seq(ORFF[p,1],ORFF[p,2])]=1
      Pdf<-rbind(Pdf,P1)
      print((colSums(Pdf[,seq(ORFF[p,1],ORFF[p,2])])))
      MCS<-max(colSums(Pdf[,seq(ORFF[p,1],ORFF[p,2])]))
      print(MCS)
      ORFFF<-(cbind(ORFF[p,],MCS))
      Pfd1<-rbind(Pfd1,ORFFF)
    }
    
    Pfd1<-Pfd1[-1,]
    inens<-c()
    for (p in seq(1,length(Pfd1[,1]))){
      inens<-c(inens,1)
    }
    
    
    for (p in seq(1,length(Pfd1[,1]))){
      color_transparent <- adjustcolor('grey', alpha.f = inens[p])
      xbottom <- Pfd1[p,1]
      ybottom <- Pfd1[p,3]+19
      xtop <- Pfd1[p,2]
      ytop<- Pfd1[p,3]+20
      xpol <- c(xbottom,xbottom,xtop,xtop)
      ypol <- c(ybottom,ytop,ytop,ybottom)
      polygon(x=cuts(t(xpol)), y=cuts(t(ypol)), col = color_transparent, border ='black')
    }
    
    
    #ADD sncRNA
    sncRNA<-read.table('out.snc.blast')
    Pfd1<-sncRNA[,c(9,10)]
    for (p in seq(1,length(Pfd1[,1]))){
      color_transparent <- adjustcolor('black', alpha.f = .5)
      xbottom <- Pfd1[p,1]
      ybottom <- 18
      xtop <- Pfd1[p,2]
      ytop<- 19
      xpol <- c(xbottom,xbottom,xtop,xtop)
      ypol <- c(ybottom,ytop,ytop,ybottom)
      polygon(x=cuts(t(xpol)), y=cuts(t(ypol)), col = color_transparent, border =NA)
    }
    
    
    
    for (p in seq(1,length(clust1cor))){
      inens<- clust1cor[p]*2
      color_transparent <- adjustcolor('orange', alpha.f = inens)
      Gm<-p
      GM<-p+1
      xbottom <- Gm
      ybottom <- 17
      xtop <- GM
      ytop<- 17.5
      xpol <- c(xbottom,xbottom,xtop,xtop)
      ypol <- c(ybottom,ytop,ytop,ybottom)
      polygon(x=cuts(t(xpol)), y=cuts(t(ypol)), col = color_transparent, border =NA)
    }
    
    for (p in seq(1,length(clust2cor))){
      inens<- clust2cor[p]*2
      color_transparent <- adjustcolor('blue', alpha.f = inens)
      Gm<-p
      GM<-p+1
      xbottom <- Gm
      ybottom <- 16.5
      xtop <- GM
      ytop<- 17
      xpol <- c(xbottom,xbottom,xtop,xtop)
      ypol <- c(ybottom,ytop,ytop,ybottom)
      polygon(x=cuts(t(xpol)), y=cuts(t(ypol)), col = color_transparent, border =NA)
    }
    
    for (p in seq(1,length(clust3cor))){
      inens<- clust3cor[p]*2
      color_transparent <- adjustcolor('green', alpha.f = inens)
      Gm<-p
      GM<-p+1
      xbottom <- Gm
      ybottom <- 16
      xtop <- GM
      ytop<- 16.5
      xpol <- c(xbottom,xbottom,xtop,xtop)
      ypol <- c(ybottom,ytop,ytop,ybottom)
      polygon(x=cuts(t(xpol)), y=cuts(t(ypol)), col = color_transparent, border =NA)
    }
    
    for (p in seq(1,length(clust4cor))){
      inens<- clust4cor[p]*2
      color_transparent <- adjustcolor('red', alpha.f = inens)
      Gm<-p
      GM<-p+1
      xbottom <- Gm
      ybottom <- 16
      xtop <- GM
      ytop<- 15.5
      xpol <- c(xbottom,xbottom,xtop,xtop)
      ypol <- c(ybottom,ytop,ytop,ybottom)
      polygon(x=cuts(t(xpol)), y=cuts(t(ypol)), col = color_transparent, border =NA)
    }
    
    #BCTV genes
    {
      color_transparent <- adjustcolor('grey0', alpha.f = 0.5)
      sncRNA<-read.table('BCTVgenes_prot.blast.out')
      Gm<-c(min(sncRNA[sncRNA[,1]=='V1',c(9,10)]))
      GM<-c(max(sncRNA[sncRNA[,1]=='V1',c(9,10)]))
      xbottom <- Gm
      ybottom <- 5
      xtop <- GM
      ytop<- 6
      xpol <- c(xbottom,xbottom,xtop,xtop)
      ypol <- c(ybottom,ytop,ytop,ybottom)
      polygon(x=cuts(t(xpol)), y=cuts(t(ypol)), col = color_transparent, border ='black')
      
      Gm<-c(min(sncRNA[sncRNA[,1]=='V2',c(9,10)]))
      GM<-c(max(sncRNA[sncRNA[,1]=='V2',c(9,10)]))
      xbottom <- Gm
      ybottom <- 6
      xtop <- GM
      ytop<- 7
      xpol <- c(xbottom,xbottom,xtop,xtop)
      ypol <- c(ybottom,ytop,ytop,ybottom)
      polygon(x=cuts(t(xpol)), y=cuts(t(ypol)), col = color_transparent, border ='black')
      
      Gm<-c(min(sncRNA[sncRNA[,1]=='V3',c(9,10)]))
      GM<-c(max(sncRNA[sncRNA[,1]=='V3',c(9,10)]))
      xbottom <- Gm
      ybottom <- 5
      xtop <- GM
      ytop<- 6
      xpol <- c(xbottom,xbottom,xtop,xtop)
      ypol <- c(ybottom,ytop,ytop,ybottom)
      polygon(x=cuts(t(xpol)), y=cuts(t(ypol)), col = color_transparent, border ='black')
      
      color_transparent <- adjustcolor('darkslategray', alpha.f = 0.5)
      Gm<-c(min(sncRNA[sncRNA[,1]=='C1',c(9,10)]))
      GM<-c(max(sncRNA[sncRNA[,1]=='C1',c(9,10)]))
      xbottom <- Gm
      ybottom <- 5
      xtop <- GM
      ytop<- 6
      xpol <- c(xbottom,xbottom,xtop,xtop)
      ypol <- c(ybottom,ytop,ytop,ybottom)
      polygon(x=cuts(t(xpol)), y=cuts(t(ypol)), col = color_transparent, border ='black')
      
      Gm<-c(min(sncRNA[sncRNA[,1]=='C2',c(9,10)]))
      GM<-c(max(sncRNA[sncRNA[,1]=='C2',c(9,10)]))
      xbottom <- Gm
      ybottom <- 6
      xtop <- GM
      ytop<- 7
      xpol <- c(xbottom,xbottom,xtop,xtop)
      ypol <- c(ybottom,ytop,ytop,ybottom)
      polygon(x=cuts(t(xpol)), y=cuts(t(ypol)), col = color_transparent, border ='black')
      
      sncRNA[sncRNA[,1]=='C3',c(9,10)]
      Gm<-c(min(sncRNA[sncRNA[,1]=='C3',c(9,10)]))
      GM<-c(max(sncRNA[sncRNA[,1]=='C3',c(9,10)]))
      xbottom <- Gm
      ybottom <- 5
      xtop <- GM
      ytop<- 6
      xpol <- c(xbottom,xbottom,xtop,xtop)
      ypol <- c(ybottom,ytop,ytop,ybottom)
      polygon(x=cuts(t(xpol)), y=cuts(t(ypol)), col = color_transparent, border ='black')
      
      sncRNA[sncRNA[,1]=='C4',c(9,10)]
      Gm<-c(min(sncRNA[sncRNA[,1]=='C4',c(9,10)]))
      GM<-c(max(sncRNA[sncRNA[,1]=='C4',c(9,10)]))
      xbottom <- Gm
      ybottom <- 6
      xtop <- GM
      ytop<- 7
      xpol <- c(xbottom,xbottom,xtop,xtop)
      ypol <- c(ybottom,ytop,ytop,ybottom)
      polygon(x=cuts(t(xpol)), y=cuts(t(ypol)), col = color_transparent, border ='black')
    }
    
    unique(DepthMD[,2])
    ROWNamova<-CLusterz[,1]
    Structure<-as.data.frame(as.vector(CLusterz[,9]))
    colnames(Structure)<-'Clusters'
    cbind(rep('C',24),Structure)
    stuct<-paste(CLusterz[,9],rep('C',19),sep = '')
    
    ##Reatin Rownames for plots later
    SAMP<-t(as.data.frame(Depth2[ROWNamova,]))
    colnames(SAMP)<-seq(1,length(SAMP[1,]))
    rownames(SAMP)<-seq(1,length(SAMP[,1]))
    SAMP1<-cbind(SAMP,rowMeans(AFV))
    SAMP2<-SAMP1[SAMP1[,20]>0,]
    SAMP3<-SAMP2[,-20]
    SAMP4<-as.matrix(SAMP)
    SAMP5<-c()
    for (i in seq(1,length(SAMP4[,1]))){
      if (length(unique(SAMP4[i,]))>1){
        SAMP5<-rbind(SAMP5,SAMP4[i,])
      }
    }
    
    NCC1t<-data.frame()
    for (i in seq(1,length(SAMP5[,1]))){
      NCC1<-cbind(rep(i,5),c(1,2,3,4,5))
      for (l in seq(1,length(SAMP5[1,]))){
        n1=0
        n2=0
        n3=0
        n4=0
        n5=0
        Nc<-c()
        SAMP7<-SAMP5[i,l]
        n0<-i
        n1<-length(SAMP7[SAMP7=='1'])
        n2<-length(SAMP7[SAMP7=='2'])
        n3<-length(SAMP7[SAMP7=='3'])
        n4<-length(SAMP7[SAMP7=='4'])
        n5<-length(SAMP7[SAMP7=='5'])
        Nc<-c(n1,n2,n3,n4,n5)
        NCC1<-cbind(NCC1,Nc)
      }
      NCC1t<-rbind(NCC1t,NCC1)
      #colnames(NCC1t)<-c('1','2','3','4','5')
    }
    
    rownames(NCC1t)<-paste(NCC1t[,1],NCC1t[,2],sep = '_')
    gsub(pattern = '_',x = rownames(NCC1t),replacement = ' ')
    NCC1t<-NCC1t[,c(-1,-2)]
    
    
    library(ade4)
    AMOVA_out<-data.frame()
    for (i in seq(1,length(NCC1t[,1]))){
      Amovdf<-NCC1t
      AMOV<-amova(samples=data.frame(Amovdf),structures = data.frame(stuct))
      print(AMOV$componentsofcovariance[1,])
      AMOVA_out<-rbind(AMOVA_out,AMOV$componentsofcovariance[1,])
    }
    
    NCC1t<-NCC1t[rowMeans(NCC1t)>0,]
    AMOV<-amova(samples=data.frame(NCC1t),structures = data.frame(stuct))
    
    ###GLM approch
    stuct<-data.frame(stuct)
    stuct<-cbind(rownames(stuct),stuct)
    colsClust<-c(1,2,3,4,5)
    NCC1tP1<-NCC1t[,colsClust]
    NCC1tP2<-NCC1t[,-colsClust]
    
    
    
    TT<-c(1,2,3,4)
    for (i in seq(1,length(NCC1tP2[,1]))){
      Tt<-t.test(NCC1tP1[i,],NCC1tP2[i,])
      Ta<-c(i,Tt$p.value,Tt$estimate)
      TT<-rbind(TT,Ta)
    }
    TTo<-TT[order(TT[,2]),]
    TTo<-TTo[TTo[,2]<0.05,]
    Ortest<-as.vector(TT[order(TT[,2]),][,1])
    NCC1t[TTo[,1],]
    TTo1<-cbind(NCC1t[TTo[,1],],TTo)
    TTo2<-TTo1[TTo1[,22]>0,]
    TTo2<-na.omit(TTo2)
    TTo1[order(TTo1$V1,decreasing = F),]
    vec_pos<-c()
    vec_id<-c()
    for (i in rownames(TTo2)){
      vec_pos<-c(vec_pos,strsplit(i,split = "_")[[1]][1])
      vec_id<-c(vec_id,strsplit(i,split = "_")[[1]][2])
    }
    
    vec_pos<-as.numeric(vec_pos)
    TTo2<-cbind(TTo2,vec_pos,vec_id)
    
    for (p in seq(1,length(TTo2[,1]))){
      color_transparent <- adjustcolor('#9F7000', alpha.f = 0.7)
      Gm<-TTo2[p,25]
      GM<-TTo2[p,25]+1
      xbottom <- Gm
      ybottom <- 16.5
      xtop <- GM
      ytop<- 17.5
      xpol <- c(xbottom,xbottom,xtop,xtop)
      ypol <- c(ybottom,ytop,ytop,ybottom)
      polygon(x=cuts(t(xpol)), y=cuts(t(ypol)), col = color_transparent, border =NA)
    }
    TTfinalTable1<-TTo2
    #cluster 2
    stuct<-data.frame(stuct)
    stuct<-cbind(rownames(stuct),stuct)
    colsClust<-c(6,7,8,9,10,11)
    NCC1tP1<-NCC1t[,colsClust]
    NCC1tP2<-NCC1t[,-colsClust]
    
    TT<-c(1,2,3,4)
    for (i in seq(1,length(NCC1tP2[,1]))){
      Tt<-t.test(NCC1tP1[i,],NCC1tP2[i,])
      Ta<-c(i,Tt$p.value,Tt$estimate)
      TT<-rbind(TT,Ta)
    }
    TTo<-TT[order(TT[,2]),]
    TTo<-TTo[TTo[,2]<0.05,]
    Ortest<-as.vector(TT[order(TT[,2]),][,1])
    NCC1t[TTo[,1],]
    TTo1<-cbind(NCC1t[TTo[,1],],TTo)
    TTo2<-TTo1[TTo1[,22]>0,]
    TTo2<-na.omit(TTo2)
    TTo1[order(TTo1$V1,decreasing = F),]
    vec_pos<-c()
    vec_id<-c()
    for (i in rownames(TTo2)){
      vec_pos<-c(vec_pos,strsplit(i,split = "_")[[1]][1])
      vec_id<-c(vec_id,strsplit(i,split = "_")[[1]][2])
    }
    
    vec_pos<-as.numeric(vec_pos)
    TTo2<-cbind(TTo2,vec_pos,vec_id)
    
    for (p in seq(1,length(TTo2[,1]))){
      color_transparent <- adjustcolor('#228B00', alpha.f = 0.7)
      Gm<-TTo2[p,25]
      GM<-TTo2[p,25]+1
      xbottom <- Gm
      ybottom <- 15.5
      xtop <- GM
      ytop<- 16.5
      xpol <- c(xbottom,xbottom,xtop,xtop)
      ypol <- c(ybottom,ytop,ytop,ybottom)
      polygon(x=cuts(t(xpol)), y=cuts(t(ypol)), col = color_transparent, border =NA)
    }
    TTfinalTable2<-TTo2
    
    #cluster 3
    stuct<-data.frame(stuct)
    stuct<-cbind(rownames(stuct),stuct)
    colsClust<-c(11,12)
    NCC1tP1<-NCC1t[,colsClust]
    NCC1tP2<-NCC1t[,-colsClust]
    
    TT<-c(1,2,3,4)
    for (i in seq(1,length(NCC1tP2[,1]))){
      Tt<-t.test(NCC1tP1[i,],NCC1tP2[i,])
      Ta<-c(i,Tt$p.value,Tt$estimate)
      TT<-rbind(TT,Ta)
    }
    TTo<-TT[order(TT[,2]),]
    TTo<-TTo[TTo[,2]<0.05,]
    Ortest<-as.vector(TT[order(TT[,2]),][,1])
    NCC1t[TTo[,1],]
    TTo1<-cbind(NCC1t[TTo[,1],],TTo)
    TTo2<-TTo1[TTo1[,22]>0,]
    TTo2<-na.omit(TTo2)
    TTo1[order(TTo1$V1,decreasing = F),]
    vec_pos<-c()
    vec_id<-c()
    for (i in rownames(TTo2)){
      vec_pos<-c(vec_pos,strsplit(i,split = "_")[[1]][1])
      vec_id<-c(vec_id,strsplit(i,split = "_")[[1]][2])
    }
    
    vec_pos<-as.numeric(vec_pos)
    TTo2<-cbind(TTo2,vec_pos,vec_id)
    
    for (p in seq(1,length(TTo2[,1]))){
      color_transparent <- adjustcolor('#009681', alpha.f = 0.7)
      Gm<-TTo2[p,25]
      GM<-TTo2[p,25]+1
      xbottom <- Gm
      ybottom <- 14.5
      xtop <- GM
      ytop<- 15.5 
      xpol <- c(xbottom,xbottom,xtop,xtop)
      ypol <- c(ybottom,ytop,ytop,ybottom)
      polygon(x=cuts(t(xpol)), y=cuts(t(ypol)), col = color_transparent, border =NA)
    }
    TTfinalTable3<-TTo2
    #cluster 4
    stuct<-data.frame(stuct)
    stuct<-cbind(rownames(stuct),stuct)
    colsClust<-c(13,15,16,17,18)
    NCC1tP1<-NCC1t[,colsClust]
    NCC1tP2<-NCC1t[,-colsClust]
    
    TT<-c(1,2,3,4)
    for (i in seq(1,length(NCC1tP2[,1]))){
      Tt<-t.test(NCC1tP1[i,],NCC1tP2[i,])
      Ta<-c(i,Tt$p.value,Tt$estimate)
      TT<-rbind(TT,Ta)
    }
    TTo<-TT[order(TT[,2]),]
    TTo<-TTo[TTo[,2]<0.05,]
    Ortest<-as.vector(TT[order(TT[,2]),][,1])
    NCC1t[TTo[,1],]
    TTo1<-cbind(NCC1t[TTo[,1],],TTo)
    TTo2<-TTo1[TTo1[,22]>0,]
    TTo2<-na.omit(TTo2)
    TTo1[order(TTo1$V1,decreasing = F),]
    vec_pos<-c()
    vec_id<-c()
    for (i in rownames(TTo2)){
      vec_pos<-c(vec_pos,strsplit(i,split = "_")[[1]][1])
      vec_id<-c(vec_id,strsplit(i,split = "_")[[1]][2])
    }
    
    vec_pos<-as.numeric(vec_pos)
    TTo2<-cbind(TTo2,vec_pos,vec_id)
    
    for (p in seq(1,length(TTo2[,1]))){
      color_transparent <- adjustcolor('#B646C7', alpha.f = 0.7)
      Gm<-TTo2[p,25]
      GM<-TTo2[p,25]+1
      xbottom <- Gm
      ybottom <- 13.5
      xtop <- GM
      ytop<- 14.5
      xpol <- c(xbottom,xbottom,xtop,xtop)
      ypol <- c(ybottom,ytop,ytop,ybottom)
      polygon(x=cuts(t(xpol)), y=cuts(t(ypol)), col = color_transparent, border =NA)
    }
    TTfinalTable4<-TTo2
    
  }       
    
  
  #Colorado
  {
    unique(DepthMD[,2])
    ROWNamova<-rownames(DepthMD[DepthMD[,2]=='BCTV-CO',])
    Structure<-as.data.frame(as.vector(DepthMD[DepthMD[,2]=='BCTV-CO',4]))
    colnames(Structure)<-'State'
    SAMP<-t(as.data.frame(Depth2[ROWNamova,]))
    colnames(SAMP)<-seq(1,length(SAMP[1,]))
    SAMP1<-cbind(SAMP,AFV[,4])
    SAMP2<-SAMP1[SAMP1[,36]>0,]
    SAMP3<-SAMP2[,-36]
    SAMP4<-as.matrix(SAMP)
    SAMP5<-c()
    for (i in seq(1,length(SAMP4[,1]))){
      if (length(unique(SAMP4[i,]))>1){
        SAMP5<-rbind(SAMP5,SAMP4[i,])
      }
    }
    
    SMAP5<-data.frame(SAMP5)
    DIST<-dist(t(SAMP5),method='euclidean')
    #install.packages('dendextend')
    library(dendextend)
    hc<-hclust(DIST)
    dhc <- as.dendrogram(hc)
    dend2 <- color_branches(dhc, k = 8)
    dend2<-set(dend2, "branches_lwd", 3)
    plot(dend2)
    get_leaves_branches_col(dend2)
    # creates a single item vector of the clusters    
    myclusters <- cutree(dhc, k=8)
    # make the dataframe of two columns cluster number and label
    clusterDF <-  data.frame(Cluster = as.numeric(unlist(myclusters)),
                             Branch = names(myclusters))
    clusterDF[order(clusterDF[,1]),]
    ROWNamova[clusterDF[order(clusterDF[,1]),2]]
    CLusterz<-cbind(DepthMD[ROWNamova[clusterDF[order(clusterDF[,1]),2]],],clusterDF[order(clusterDF[,1]),])
    unique(CLusterz[,9])
    
    
    FstALL<-data.frame()
    Cor1<-data.frame()
    AFV<-data.frame(1,2,3,4)
    colnames(AFV)<-c('1','2','3','4')
    clust1<-CLusterz[CLusterz[,9]=='6',]
    clust2<-CLusterz[CLusterz[,9]=='3',]
    clust3<-CLusterz[CLusterz[,9]=='8',]
    clust4<-CLusterz[CLusterz[,9]=='9',]
    clust5<-CLusterz[CLusterz[,9]=='7',]
    clust6<-CLusterz[CLusterz[,9]=='6',]
    
    for (i in seq(1,3175)){
      i<-i
      print(i)
      X<-as.vector(Depth2[rownames(clust1),i])
      Xa<-(length(X[X==1])/length(X))^2
      Xt<-(length(X[X==2])/length(X))^2
      Xg<-(length(X[X==3])/length(X))^2
      Xc<-(length(X[X==4])/length(X))^2
      Xx<-(length(X[X==5])/length(X))^2
      SX<-1-sum(Xa,Xt,Xg,Xc,Xx)
      V1<-c(Xa,Xt,Xg,Xc,Xx)
      SUMX<-(SX)
      
      X1<-as.vector(Depth2[rownames(clust2),i])
      Xa1<-(length(X1[X1==1])/length(X1))^2
      Xt1<-(length(X1[X1==2])/length(X1))^2
      Xg1<-(length(X1[X1==3])/length(X1))^2
      Xc1<-(length(X1[X1==4])/length(X1))^2
      Xx1<-(length(X1[X1==5])/length(X1))^2
      SX1<-1-sum(Xa1,Xt1,Xg1,Xc1,Xx1)
      SUMX1<-(SX1)
      FstSite<-((SX+SX1)-(SX1))/(SX+SX1)
      V2<-c(Xa1,Xt1,Xg1,Xc1,Xx1)
      C1<-cor(V1,V2)
      
      X2<-as.vector(Depth2[rownames(clust3),i])
      Xa2<-(length(X2[X2==1])/length(X2))^2
      Xt2<-(length(X2[X2==2])/length(X2))^2
      Xg2<-(length(X2[X2==3])/length(X2))^2
      Xc2<-(length(X2[X2==4])/length(X2))^2
      Xx2<-(length(X2[X2==5])/length(X2))^2
      SX2<-1-sum(Xa2,Xt2,Xg2,Xc2,Xx2)
      V3<-c(Xa2,Xt2,Xg2,Xc2,Xx2)
      
      X3<-as.vector(Depth2[rownames(clust4),i])
      Xa3<-(length(X3[X3==1])/length(X3))^2
      Xt3<-(length(X3[X3==2])/length(X3))^2
      Xg3<-(length(X3[X3==3])/length(X3))^2
      Xc3<-(length(X3[X3==4])/length(X3))^2
      Xx3<-(length(X3[X3==5])/length(X3))^2
      SX3<-1-sum(Xa3,Xt3,Xg3,Xc3,Xx3)
      V4<-c(Xa3,Xt3,Xg3,Xc3,Xx3)
      
      
      C12id<-cor(V1,V2)
      C13id<-cor(V1,V3)
      C14id<-cor(V1,V4)
      
      
      C21mt<-cor(V2,V1)
      C23mt<-cor(V2,V3)
      C24mt<-cor(V2,V4)
      
      
      C31or<-cor(V3,V1)
      C32or<-cor(V3,V2)
      C34or<-cor(V3,V4)
      
      
      C41wy<-cor(V4,V1)
      C42wy<-cor(V4,V2)
      C43wy<-cor(V4,V3)
      
      
      C1cor<-c(C12id,C13id,C14id,C21mt,C23mt,C24mt,
               C31or,C32or,C34or,C41wy,C42wy,C43wy)
      
      AFV1<-c(SX,SX1,SX2,SX3)
      AFV<-rbind(AFV,AFV1)
      Cor1<-rbind(Cor1,C1cor)
    }
    AFV<-AFV[-1,]
    colMeans(AFV)
    
    Cor1<-(1-abs(Cor1))
    
    
    Cor1
    clust1cor<-rowMeans(Cor1[,c(1,2,3)])
    clust2cor<-rowMeans(Cor1[,c(4,5,6)])
    clust3cor<-rowMeans(Cor1[,c(7,8,9)])
    clust4cor<-rowMeans(Cor1[,c(10,11,12)])
    
    clust1cor[clust1cor<0.2]=0
    clust2cor[clust2cor<0.2]=0
    clust3cor[clust3cor<0.2]=0
    clust4cor[clust4cor<0.2]=0
    
    colMeans(cbind(clust1cor,clust2cor,clust3cor,clust4cor))
    
    
    par(mfrow=c(1,1))
    plot(NULL,ylim=c(0,40),xlim = c(0,3175))
    
    ORFF<-ORFs[ORFs[,1]=='giCO',c(13,14)]
    
    lorf<-c()
    for (i in seq(1,length(ORFF[,1]))){
      print(max(c(ORFF[i,1],ORFF[i,2]))-min(c(ORFF[i,1],ORFF[i,2])))
      lorf<-c(lorf,max(c(ORFF[i,1],ORFF[i,2]))-min(c(ORFF[i,1],ORFF[i,2])))
    }
    olorf<-order(lorf)
    oolorf<-olorf[seq(1,length(lorf)-7)]
    ORFF<-ORFF[oolorf,]
    
    
    lorf<-c()
    for (p in seq(1,length(ORFF[,1]))){
      ll<-length(seq(ORFF[p,1],ORFF[p,2]))
      lorf<-c(lorf,ll)
    }
    
    lorf<-c()
    for (p in seq(1,length(ORFF[,1]))){
      ll<-(length(seq(ORFF[p,1],ORFF[p,2])))
      lorf<-c(lorf,ll)
      
    }
    lorf
    lorfo<-order(lorf,decreasing = F)
    
    Pdf<-c(rep(0,3175))
    Pfd1<-c(0,0,0)
    for (p in lorfo){
      print(ORFF[p,])
      P1<-c(rep(0,3175))
      P1[seq(ORFF[p,1],ORFF[p,2])]=1
      Pdf<-rbind(Pdf,P1)
      print((colSums(Pdf[,seq(ORFF[p,1],ORFF[p,2])])))
      MCS<-max(colSums(Pdf[,seq(ORFF[p,1],ORFF[p,2])]))
      print(MCS)
      ORFFF<-(cbind(ORFF[p,],MCS))
      Pfd1<-rbind(Pfd1,ORFFF)
    }
    
    Pfd1<-Pfd1[-1,]
    inens<-c()
    for (p in seq(1,length(Pfd1[,1]))){
      inens<-c(inens,1)
    }
    
    
    for (p in seq(1,length(Pfd1[,1]))){
      color_transparent <- adjustcolor('grey', alpha.f = inens[p])
      xbottom <- Pfd1[p,1]
      ybottom <- Pfd1[p,3]+19
      xtop <- Pfd1[p,2]
      ytop<- Pfd1[p,3]+20
      xpol <- c(xbottom,xbottom,xtop,xtop)
      ypol <- c(ybottom,ytop,ytop,ybottom)
      polygon(x=cuts(t(xpol)), y=cuts(t(ypol)), col = color_transparent, border ='black')
    }
    
    
    #ADD sncRNA
    sncRNA<-read.table('out.snc.blast')
    Pfd1<-sncRNA[,c(9,10)]
    for (p in seq(1,length(Pfd1[,1]))){
      color_transparent <- adjustcolor('black', alpha.f = .5)
      xbottom <- Pfd1[p,1]
      ybottom <- 18
      xtop <- Pfd1[p,2]
      ytop<- 19
      xpol <- c(xbottom,xbottom,xtop,xtop)
      ypol <- c(ybottom,ytop,ytop,ybottom)
      polygon(x=cuts(t(xpol)), y=cuts(t(ypol)), col = color_transparent, border =NA)
    }
    
    
    
    for (p in seq(1,length(clust1cor))){
      inens<- clust1cor[p]*2
      color_transparent <- adjustcolor('orange', alpha.f = inens)
      Gm<-p
      GM<-p+1
      xbottom <- Gm
      ybottom <- 17
      xtop <- GM
      ytop<- 17.5
      xpol <- c(xbottom,xbottom,xtop,xtop)
      ypol <- c(ybottom,ytop,ytop,ybottom)
      polygon(x=cuts(t(xpol)), y=cuts(t(ypol)), col = color_transparent, border =NA)
    }
    
    for (p in seq(1,length(clust2cor))){
      inens<- clust2cor[p]*2
      color_transparent <- adjustcolor('blue', alpha.f = inens)
      Gm<-p
      GM<-p+1
      xbottom <- Gm
      ybottom <- 16.5
      xtop <- GM
      ytop<- 17
      xpol <- c(xbottom,xbottom,xtop,xtop)
      ypol <- c(ybottom,ytop,ytop,ybottom)
      polygon(x=cuts(t(xpol)), y=cuts(t(ypol)), col = color_transparent, border =NA)
    }
    
    for (p in seq(1,length(clust3cor))){
      inens<- clust3cor[p]*2
      color_transparent <- adjustcolor('green', alpha.f = inens)
      Gm<-p
      GM<-p+1
      xbottom <- Gm
      ybottom <- 16
      xtop <- GM
      ytop<- 16.5
      xpol <- c(xbottom,xbottom,xtop,xtop)
      ypol <- c(ybottom,ytop,ytop,ybottom)
      polygon(x=cuts(t(xpol)), y=cuts(t(ypol)), col = color_transparent, border =NA)
    }
    
    for (p in seq(1,length(clust4cor))){
      inens<- clust4cor[p]*2
      color_transparent <- adjustcolor('red', alpha.f = inens)
      Gm<-p
      GM<-p+1
      xbottom <- Gm
      ybottom <- 16
      xtop <- GM
      ytop<- 15.5
      xpol <- c(xbottom,xbottom,xtop,xtop)
      ypol <- c(ybottom,ytop,ytop,ybottom)
      polygon(x=cuts(t(xpol)), y=cuts(t(ypol)), col = color_transparent, border =NA)
    }
    
    #BCTV genes
    {
      color_transparent <- adjustcolor('grey0', alpha.f = 0.5)
      sncRNA<-read.table('BCTVgenes_prot.blast.out')
      Gm<-c(min(sncRNA[sncRNA[,1]=='V1',c(9,10)]))
      GM<-c(max(sncRNA[sncRNA[,1]=='V1',c(9,10)]))
      xbottom <- Gm
      ybottom <- 5
      xtop <- GM
      ytop<- 6
      xpol <- c(xbottom,xbottom,xtop,xtop)
      ypol <- c(ybottom,ytop,ytop,ybottom)
      polygon(x=cuts(t(xpol)), y=cuts(t(ypol)), col = color_transparent, border ='black')
      
      Gm<-c(min(sncRNA[sncRNA[,1]=='V2',c(9,10)]))
      GM<-c(max(sncRNA[sncRNA[,1]=='V2',c(9,10)]))
      xbottom <- Gm
      ybottom <- 6
      xtop <- GM
      ytop<- 7
      xpol <- c(xbottom,xbottom,xtop,xtop)
      ypol <- c(ybottom,ytop,ytop,ybottom)
      polygon(x=cuts(t(xpol)), y=cuts(t(ypol)), col = color_transparent, border ='black')
      
      Gm<-c(min(sncRNA[sncRNA[,1]=='V3',c(9,10)]))
      GM<-c(max(sncRNA[sncRNA[,1]=='V3',c(9,10)]))
      xbottom <- Gm
      ybottom <- 5
      xtop <- GM
      ytop<- 6
      xpol <- c(xbottom,xbottom,xtop,xtop)
      ypol <- c(ybottom,ytop,ytop,ybottom)
      polygon(x=cuts(t(xpol)), y=cuts(t(ypol)), col = color_transparent, border ='black')
      
      color_transparent <- adjustcolor('darkslategray', alpha.f = 0.5)
      Gm<-c(min(sncRNA[sncRNA[,1]=='C1',c(9,10)]))
      GM<-c(max(sncRNA[sncRNA[,1]=='C1',c(9,10)]))
      xbottom <- Gm
      ybottom <- 5
      xtop <- GM
      ytop<- 6
      xpol <- c(xbottom,xbottom,xtop,xtop)
      ypol <- c(ybottom,ytop,ytop,ybottom)
      polygon(x=cuts(t(xpol)), y=cuts(t(ypol)), col = color_transparent, border ='black')
      
      Gm<-c(min(sncRNA[sncRNA[,1]=='C2',c(9,10)]))
      GM<-c(max(sncRNA[sncRNA[,1]=='C2',c(9,10)]))
      xbottom <- Gm
      ybottom <- 6
      xtop <- GM
      ytop<- 7
      xpol <- c(xbottom,xbottom,xtop,xtop)
      ypol <- c(ybottom,ytop,ytop,ybottom)
      polygon(x=cuts(t(xpol)), y=cuts(t(ypol)), col = color_transparent, border ='black')
      
      sncRNA[sncRNA[,1]=='C3',c(9,10)]
      Gm<-c(min(sncRNA[sncRNA[,1]=='C3',c(9,10)]))
      GM<-c(max(sncRNA[sncRNA[,1]=='C3',c(9,10)]))
      xbottom <- Gm
      ybottom <- 5
      xtop <- GM
      ytop<- 6
      xpol <- c(xbottom,xbottom,xtop,xtop)
      ypol <- c(ybottom,ytop,ytop,ybottom)
      polygon(x=cuts(t(xpol)), y=cuts(t(ypol)), col = color_transparent, border ='black')
      
      sncRNA[sncRNA[,1]=='C4',c(9,10)]
      Gm<-c(min(sncRNA[sncRNA[,1]=='C4',c(9,10)]))
      GM<-c(max(sncRNA[sncRNA[,1]=='C4',c(9,10)]))
      xbottom <- Gm
      ybottom <- 6
      xtop <- GM
      ytop<- 7
      xpol <- c(xbottom,xbottom,xtop,xtop)
      ypol <- c(ybottom,ytop,ytop,ybottom)
      polygon(x=cuts(t(xpol)), y=cuts(t(ypol)), col = color_transparent, border ='black')
    }
    
    unique(DepthMD[,2])
    ROWNamova<-CLusterz[,1]
    Structure<-as.data.frame(as.vector(CLusterz[,9]))
    colnames(Structure)<-'Clusters'
    cbind(rep('C',19),Structure)
    stuct<-paste(CLusterz[,9],rep('C',19),sep = '')
    
    ##Reatin Rownames for plots later
    SAMP<-t(as.data.frame(Depth2[ROWNamova,]))
    colnames(SAMP)<-seq(1,length(SAMP[1,]))
    rownames(SAMP)<-seq(1,length(SAMP[,1]))
    SAMP1<-cbind(SAMP,rowMeans(AFV))
    SAMP2<-SAMP1[SAMP1[,20]>0,]
    SAMP3<-SAMP2[,-20]
    SAMP4<-as.matrix(SAMP)
    SAMP5<-c()
    for (i in seq(1,length(SAMP4[,1]))){
      if (length(unique(SAMP4[i,]))>1){
        SAMP5<-rbind(SAMP5,SAMP4[i,])
      }
    }
    
    NCC1t<-data.frame()
    for (i in seq(1,length(SAMP5[,1]))){
      NCC1<-cbind(rep(i,5),c(1,2,3,4,5))
      for (l in seq(1,length(SAMP5[1,]))){
        n1=0
        n2=0
        n3=0
        n4=0
        n5=0
        Nc<-c()
        SAMP7<-SAMP5[i,l]
        n0<-i
        n1<-length(SAMP7[SAMP7=='1'])
        n2<-length(SAMP7[SAMP7=='2'])
        n3<-length(SAMP7[SAMP7=='3'])
        n4<-length(SAMP7[SAMP7=='4'])
        n5<-length(SAMP7[SAMP7=='5'])
        Nc<-c(n1,n2,n3,n4,n5)
        NCC1<-cbind(NCC1,Nc)
      }
      NCC1t<-rbind(NCC1t,NCC1)
      #colnames(NCC1t)<-c('1','2','3','4','5')
    }
    
    rownames(NCC1t)<-paste(NCC1t[,1],NCC1t[,2],sep = '_')
    gsub(pattern = '_',x = rownames(NCC1t),replacement = ' ')
    NCC1t<-NCC1t[,c(-1,-2)]
    
    
    library(ade4)
    AMOVA_out<-data.frame()
    for (i in seq(1,length(NCC1t[,1]))){
      Amovdf<-NCC1t
      AMOV<-amova(samples=data.frame(Amovdf),structures = data.frame(stuct))
      print(AMOV$componentsofcovariance[1,])
      AMOVA_out<-rbind(AMOVA_out,AMOV$componentsofcovariance[1,])
    }
    
    NCC1t<-NCC1t[rowMeans(NCC1t)>0,]
    AMOV<-amova(samples=data.frame(NCC1t),structures = data.frame(stuct))
    
    ###GLM approch
    stuct<-data.frame(stuct)
    stuct<-cbind(rownames(stuct),stuct)
    colsClust<-c(3,4,5,6,10,11,12,13,14,15,16,17,18,19,20,21)
    NCC1tP1<-NCC1t[,colsClust]
    NCC1tP2<-NCC1t[,-colsClust]
    
    
    
    TT<-c(1,2,3,4)
    for (i in seq(1,length(NCC1tP2[,1]))){
      Tt<-t.test(NCC1tP1[i,],NCC1tP2[i,])
      Ta<-c(i,Tt$p.value,Tt$estimate)
      TT<-rbind(TT,Ta)
    }
    TTo<-TT[order(TT[,2]),]
    TTo<-TTo[TTo[,2]<0.05,]
    Ortest<-as.vector(TT[order(TT[,2]),][,1])
    NCC1t[TTo[,1],]
    TTo1<-cbind(NCC1t[TTo[,1],],TTo)
    TTo2<-TTo1[TTo1[,22]>0,]
    TTo2<-na.omit(TTo2)
    TTo1[order(TTo1$V1,decreasing = F),]
    vec_pos<-c()
    vec_id<-c()
    for (i in rownames(TTo2)){
      vec_pos<-c(vec_pos,strsplit(i,split = "_")[[1]][1])
      vec_id<-c(vec_id,strsplit(i,split = "_")[[1]][2])
    }
    
    vec_pos<-as.numeric(vec_pos)
    TTo2<-cbind(TTo2,vec_pos,vec_id)
    
    for (p in seq(1,length(TTo2[,1]))){
      color_transparent <- adjustcolor('#767F00', alpha.f = 0.7)
      Gm<-TTo2[p,36]
      GM<-TTo2[p,36]+1
      xbottom <- Gm
      ybottom <- 16.5
      xtop <- GM
      ytop<- 17.5
      xpol <- c(xbottom,xbottom,xtop,xtop)
      ypol <- c(ybottom,ytop,ytop,ybottom)
      polygon(x=cuts(t(xpol)), y=cuts(t(ypol)), col = color_transparent, border =NA)
    }
    TTfinalTable1<-TTo2
    #cluster 2
    stuct<-data.frame(stuct)
    stuct<-cbind(rownames(stuct),stuct)
    colsClust<-c(23,28,29,30,12)
    NCC1tP1<-NCC1t[,colsClust]
    NCC1tP2<-NCC1t[,-colsClust]
    
    TT<-c(1,2,3,4)
    for (i in seq(1,length(NCC1tP2[,1]))){
      Tt<-t.test(NCC1tP1[i,],NCC1tP2[i,])
      Ta<-c(i,Tt$p.value,Tt$estimate)
      TT<-rbind(TT,Ta)
    }
    TTo<-TT[order(TT[,2]),]
    TTo<-TTo[TTo[,2]<0.05,]
    Ortest<-as.vector(TT[order(TT[,2]),][,1])
    NCC1t[TTo[,1],]
    TTo1<-cbind(NCC1t[TTo[,1],],TTo)
    TTo2<-TTo1[TTo1[,22]>0,]
    TTo2<-na.omit(TTo2)
    TTo1[order(TTo1$V1,decreasing = F),]
    vec_pos<-c()
    vec_id<-c()
    for (i in rownames(TTo2)){
      vec_pos<-c(vec_pos,strsplit(i,split = "_")[[1]][1])
      vec_id<-c(vec_id,strsplit(i,split = "_")[[1]][2])
    }
    
    vec_pos<-as.numeric(vec_pos)
    TTo2<-cbind(TTo2,vec_pos,vec_id)
    
    for (p in seq(1,length(TTo2[,1]))){
      color_transparent <- adjustcolor('#008F00', alpha.f = 0.7)
      Gm<-TTo2[p,36]
      GM<-TTo2[p,36]+1
      xbottom <- Gm
      ybottom <- 15.5
      xtop <- GM
      ytop<- 16.5
      xpol <- c(xbottom,xbottom,xtop,xtop)
      ypol <- c(ybottom,ytop,ytop,ybottom)
      polygon(x=cuts(t(xpol)), y=cuts(t(ypol)), col = color_transparent, border =NA)
    }
    TTfinalTable2<-TTo2
    
    #cluster 3
    stuct<-data.frame(stuct)
    stuct<-cbind(rownames(stuct),stuct)
    colsClust<-c(25,26,32,34)
    NCC1tP1<-NCC1t[,colsClust]
    NCC1tP2<-NCC1t[,-colsClust]
    
    TT<-c(1,2,3,4)
    for (i in seq(1,length(NCC1tP2[,1]))){
      Tt<-t.test(NCC1tP1[i,],NCC1tP2[i,])
      Ta<-c(i,Tt$p.value,Tt$estimate)
      TT<-rbind(TT,Ta)
    }
    TTo<-TT[order(TT[,2]),]
    TTo<-TTo[TTo[,2]<0.05,]
    Ortest<-as.vector(TT[order(TT[,2]),][,1])
    NCC1t[TTo[,1],]
    TTo1<-cbind(NCC1t[TTo[,1],],TTo)
    TTo2<-TTo1[TTo1[,22]>0,]
    TTo2<-na.omit(TTo2)
    TTo1[order(TTo1$V1,decreasing = F),]
    vec_pos<-c()
    vec_id<-c()
    for (i in rownames(TTo2)){
      vec_pos<-c(vec_pos,strsplit(i,split = "_")[[1]][1])
      vec_id<-c(vec_id,strsplit(i,split = "_")[[1]][2])
    }
    
    vec_pos<-as.numeric(vec_pos)
    TTo2<-cbind(TTo2,vec_pos,vec_id)
    
    for (p in seq(1,length(TTo2[,1]))){
      color_transparent <- adjustcolor('#009681', alpha.f = 0.7)
      Gm<-TTo2[p,36]
      GM<-TTo2[p,36]+1
      xbottom <- Gm
      ybottom <- 14.5
      xtop <- GM
      ytop<- 15.5 
      xpol <- c(xbottom,xbottom,xtop,xtop)
      ypol <- c(ybottom,ytop,ytop,ybottom)
      polygon(x=cuts(t(xpol)), y=cuts(t(ypol)), col = color_transparent, border =NA)
    }
    TTfinalTable3<-TTo2
    #cluster 4
    stuct<-data.frame(stuct)
    stuct<-cbind(rownames(stuct),stuct)
    colsClust<-c(22,24,27,35)
    NCC1tP1<-NCC1t[,colsClust]
    NCC1tP2<-NCC1t[,-colsClust]
    
    TT<-c(1,2,3,4)
    for (i in seq(1,length(NCC1tP2[,1]))){
      Tt<-t.test(NCC1tP1[i,],NCC1tP2[i,])
      Ta<-c(i,Tt$p.value,Tt$estimate)
      TT<-rbind(TT,Ta)
    }
    TTo<-TT[order(TT[,2]),]
    TTo<-TTo[TTo[,2]<0.05,]
    Ortest<-as.vector(TT[order(TT[,2]),][,1])
    NCC1t[TTo[,1],]
    TTo1<-cbind(NCC1t[TTo[,1],],TTo)
    TTo2<-TTo1[TTo1[,22]>0,]
    TTo2<-na.omit(TTo2)
    TTo1[order(TTo1$V1,decreasing = F),]
    vec_pos<-c()
    vec_id<-c()
    for (i in rownames(TTo2)){
      vec_pos<-c(vec_pos,strsplit(i,split = "_")[[1]][1])
      vec_id<-c(vec_id,strsplit(i,split = "_")[[1]][2])
    }
    
    vec_pos<-as.numeric(vec_pos)
    TTo2<-cbind(TTo2,vec_pos,vec_id)
    
    for (p in seq(1,length(TTo2[,1]))){
      color_transparent <- adjustcolor('#008BC1', alpha.f = 0.7)
      Gm<-TTo2[p,36]
      GM<-TTo2[p,36]+1
      xbottom <- Gm
      ybottom <- 13.5
      xtop <- GM
      ytop<- 14.5
      xpol <- c(xbottom,xbottom,xtop,xtop)
      ypol <- c(ybottom,ytop,ytop,ybottom)
      polygon(x=cuts(t(xpol)), y=cuts(t(ypol)), col = color_transparent, border =NA)
    }
    TTfinalTable4<-TTo2
    
  }    
    
    
    
    
    
    
    
    
    
    
    
    
    AMOVA_out<-data.frame()
    for (i in seq(1,length(NCC1t[,1]))){
      Amovdf<-NCC1t
      AMOV<-amova(samples=data.frame(Amovdf),structures = data.frame(stuct))
      print(AMOV$componentsofcovariance[1,])
      AMOVA_out<-rbind(AMOVA_out,AMOV$componentsofcovariance[1,])
      
    }
    NCC1t<-NCC1t[rowMeans(NCC1t)>0,]
    AMOV<-amova(samples=data.frame(NCC1t),structures = data.frame(stuct))
    NCC1t
    
    NCC1<-data.frame()
    for (i in seq(1,length(SAMP5[,1]))){
      n1=0
      n2=0
      n3=0
      n4=0
      n5=0
      Nc<-c()
      SAMP6<-SAMP5[,CLusterz[CLusterz[,9]=='1',10]]
      SAMP7<-SAMP6[i,]
      n1<-length(SAMP7[SAMP7=='1'])
      n2<-length(SAMP7[SAMP7=='2'])
      n3<-length(SAMP7[SAMP7=='3'])
      n4<-length(SAMP7[SAMP7=='4'])
      n5<-length(SAMP7[SAMP7=='5'])
      Nc<-c(n1,n2,n3,n4,n5)
      NCC1<-rbind(NCC1,Nc)
      colnames(NCC1)<-c('1','2','3','4','5')
    }
    
    NCC2<-data.frame()
    for (i in seq(1,length(SAMP5[,1]))){
      n1=0
      n2=0
      n3=0
      n4=0
      n5=0
      Nc<-c()
      SAMP6<-SAMP5[,CLusterz[CLusterz[,9]=='2',10]]
      SAMP7<-SAMP6[i,]
      n1<-length(SAMP7[SAMP7=='1'])
      n2<-length(SAMP7[SAMP7=='2'])
      n3<-length(SAMP7[SAMP7=='3'])
      n4<-length(SAMP7[SAMP7=='4'])
      n5<-length(SAMP7[SAMP7=='5'])
      Nc<-c(n1,n2,n3,n4,n5)
      NCC2<-rbind(NCC2,Nc)
      colnames(NCC2)<-c('1','2','3','4','5')
    }
    NCC3<-data.frame()
    for (i in seq(1,length(SAMP5[,1]))){
      n1=0
      n2=0
      n3=0
      n4=0
      n5=0
      Nc<-c()
      SAMP6<-SAMP5[,CLusterz[CLusterz[,9]=='3',10]]
      SAMP7<-SAMP6[i,]
      n1<-length(SAMP7[SAMP7=='1'])
      n2<-length(SAMP7[SAMP7=='2'])
      n3<-length(SAMP7[SAMP7=='3'])
      n4<-length(SAMP7[SAMP7=='4'])
      n5<-length(SAMP7[SAMP7=='5'])
      Nc<-c(n1,n2,n3,n4,n5)
      NCC3<-rbind(NCC3,Nc)
      colnames(NCC3)<-c('1','2','3','4','5')
    }
    
    NCC4<-data.frame()
    for (i in seq(1,length(SAMP5[,1]))){
      n1=0
      n2=0
      n3=0
      n4=0
      n5=0
      Nc<-c()
      SAMP6<-SAMP5[,CLusterz[CLusterz[,9]=='4',10]]
      SAMP7<-SAMP6[i,]
      n1<-length(SAMP7[SAMP7=='1'])
      n2<-length(SAMP7[SAMP7=='2'])
      n3<-length(SAMP7[SAMP7=='3'])
      n4<-length(SAMP7[SAMP7=='4'])
      n5<-length(SAMP7[SAMP7=='5'])
      Nc<-c(n1,n2,n3,n4,n5)
      NCC4<-rbind(NCC4,Nc)
      colnames(NCC4)<-c('1','2','3','4','5')
    }
    
    NCC5<-data.frame()
    for (i in seq(1,length(SAMP5[,1]))){
      n1=0
      n2=0
      n3=0
      n4=0
      n5=0
      Nc<-c()
      SAMP6<-SAMP5[,CLusterz[CLusterz[,9]=='5',10]]
      SAMP7<-SAMP6[i,]
      n1<-length(SAMP7[SAMP7=='1'])
      n2<-length(SAMP7[SAMP7=='2'])
      n3<-length(SAMP7[SAMP7=='3'])
      n4<-length(SAMP7[SAMP7=='4'])
      n5<-length(SAMP7[SAMP7=='5'])
      Nc<-c(n1,n2,n3,n4,n5)
      NCC5<-rbind(NCC5,Nc)
      colnames(NCC5)<-c('1','2','3','4','5')
    }
    
    NCC6<-data.frame()
    for (i in seq(1,length(SAMP5[,1]))){
      n1=0
      n2=0
      n3=0
      n4=0
      n5=0
      Nc<-c()
      SAMP6<-SAMP5[,CLusterz[CLusterz[,9]=='6',10]]
      SAMP7<-SAMP6[i,]
      n1<-length(SAMP7[SAMP7=='1'])
      n2<-length(SAMP7[SAMP7=='2'])
      n3<-length(SAMP7[SAMP7=='3'])
      n4<-length(SAMP7[SAMP7=='4'])
      n5<-length(SAMP7[SAMP7=='5'])
      Nc<-c(n1,n2,n3,n4,n5)
      NCC6<-rbind(NCC6,Nc)
      colnames(NCC6)<-c('1','2','3','4','5')
    }
    
    
    
    
    install.packages('ade4')
    library(ade4)
    stuct<-c('c1','c2','c3','c5')
    AMOVA_out<-data.frame()
    for (i in seq(1,length(NCC1[,1]))){
      Amovdf<-cbind(t(NCC1[i,]),t(NCC2[i,]),t(NCC3[i,]),t(NCC5[i,]))
      colnames(Amovdf)<-c('c1','c2','c3','c5')
      AMOV<-amova(samples=data.frame(Amovdf),structures = data.frame(stuct))
      AMOVA_out<-rbind(AMOVA_out,AMOV$componentsofcovariance[1,])
      
    }
    
    
    cbind(t(NCC1[1,]),t(NCC2[1,]),t(NCC3[1,]),t(NCC5[1,]))
    
    SAMP5<-data.frame(SAMP5)
    AMOV<-amova(samples=data.frame(SAMP5),structures = data.frame(stuct))
    AMOV$results
    AMOV$componentsofcovariance
    AMOV$statphi
    
    amova()
    
    
  }
  
  
