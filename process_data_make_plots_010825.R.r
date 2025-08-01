######################################################################################################################################################
# 
# Program to make plots for original submitted manuscript. 
# 
# 01/08/2025
# 
######################################################################################################################################################



library(ncdf4)
library(fitdistrplus)
library(evd)
library(rgdal)
library(maps)
library(maptools)
library(data.table)
library(dplR)


ipdir.ip<- "D:/***/"

dir.interim<- "E:/***/"
ipdir.hist<- paste0(dir.interim, "HIST_DAMIP/sfcewind/")
ipdir.ghg<- paste0(dir.interim, "HIST_GHG/sfcewind/")
ipdir.nat<- paste0(dir.interim, "HIST_Nat/sfcewind/")
ipdir.aero<- paste0(dir.interim, "HIST_Aer/sfcewind/")
ipdir.ctl<- paste0(dir.interim, "piControl/sfcewind/")

opdir.plots<- "D:/***/"



###  1. Read in setup information such as CMIP6 model names and info

ip.fn<- paste0(ipdir.ip,'CMIP6_model_names_DAMIP_windmax_v3.csv')
ip.data<- read.table(ip.fn, header=T, sep=',')
model.id.all<- as.character(ip.data$model); hist.runid.all<- as.character(ip.data$Hist_run); hist.grid.all<- as.character(ip.data$Hist_grid)
ghg.runid.all<- as.character(ip.data$Hist_GHG_run); ghg.grid.all<- as.character(ip.data$Hist_GHG_grid)
nat.runid.all<- as.character(ip.data$Hist_Nat_run); nat.grid.all<- as.character(ip.data$Hist_Nat_grid)
aero.runid.all<- as.character(ip.data$Hist_Aer_run); aero.grid.all<- as.character(ip.data$Hist_Aer_grid)

num.cmip6.models<- length(model.id.all)

ip.fn<- paste0(ipdir.ip,'CMIP6_model_names_DAMIP_windmax_CTL_v3.csv')
ip.data<- read.table(ip.fn, header=T, sep=',')
model.id.ctl<- as.character(ip.data$model); ctl.runid.all<- as.character(ip.data$CTL_run); ctl.grid.all<- as.character(ip.data$CTL_grid)
num.ctl.models<- length(model.id.ctl)

num.yrs.ctl.model<- c(500, 600, 500, 500, 600, 250)



# read in grid information

ip.fn<- paste0(ipdir.ip,'universal_grid_EU_0p5_v2.csv')
ip.data<- read.table(ip.fn, header=T, sep=',')
grid.lat<- ip.data$LatC;  grid.lon<- ip.data$LonC;  grid.pop<- ip.data$pop;  grid.ctry<- as.character(ip.data$country)
ptr.c<- which(grid.ctry == 'PL' | grid.ctry == 'CZ' | grid.ctry == 'SK'); if (length(ptr.c) > 0) { grid.pop[ptr.c]<- ip.data$pop[ptr.c] * 0.5 }  #  adjust for low insurance uptake in eastern EU
num.grid.pts<- length(grid.lat)
all.ctry<- c(unique(grid.ctry), 'EU'); num.ctry<- length(all.ctry)



### set filtering parameters ('filt.band' contains periods being passed, 'usr.filt.n' is the order of the Butterworth filter)
filt.period<- 20
usr.filt.n<- 2





######################################################################################################################################################


###  read in observed losses and NAO, to make Figure 1


filt.period.init<- 10
usr.filt.n<- 2

ev.loss.fn<- 'D:/***/**.csv'
ip.dat<- read.table(ev.loss.fn, header=T, sep=',')
px<- which(as.integer(ip.dat$Year) >= 1950)
ip.date<- as.integer(ip.dat$Year[px])
euws.ann.loss<- as.numeric(ip.dat$EU[px])*1e-6

euws.ann.loss.anom<- euws.ann.loss - mean(euws.ann.loss)
euws.loss.10yr<- pass.filt(euws.ann.loss.anom, filt.period.init, type="low", method="Butterworth", n=usr.filt.n)


###
euws.loss.20yr<- pass.filt(euws.ann.loss, filt.period, type="low", method="Butterworth", n=usr.filt.n)
temparr1<- (euws.ann.loss - euws.loss.20yr)/euws.loss.20yr
temparr2<- (euws.ann.loss - mean(euws.ann.loss))/mean(euws.ann.loss)
sd(temparr1) / sqrt(20)

###

nao.fn<- 'D:/***/NAO_monthly_stations_UEA.txt'
ip.dat<- read.table(nao.fn, header=F, sep='	')
ip.yr<- as.integer(ip.dat$V1)

ip.NAO.ann.init<- array(0.0, 203)
for (y in 1:203) {
  temparr<- c( as.numeric(ip.dat[y, 2:4]), as.numeric(ip.dat[y, 11:13]) )
  ip.NAO.ann.init[y]<- mean(temparr)
}

ip.NAO.ann<- ip.NAO.ann.init[130:203]  #  1950-2023

NAO.10yr<- pass.filt(ip.NAO.ann, filt.period.init, type="low", method="Butterworth", n=usr.filt.n)

euws.loss.10yr.norm<- (euws.loss.10yr-mean(euws.loss.10yr)) / sd(euws.loss.10yr)
NAO.10yr.norm<- (NAO.10yr-mean(NAO.10yr)) / sd(NAO.10yr)



yrange=c(-2.5,2.544); xrange=c(1951.2,2021.7)
fn<- paste(opdir.plots,'Figure_1.png',sep='')
png(file = fn, width = 2250., height = 1000., type = 'cairo')
op <- par( mai = c(1.0, 2.05, 0.4, 0.5), cex=2.0); xx<- 1950:2023
  plot(xx, euws.loss.10yr.norm, type="l", lty = 1, cex=1.5, axes = FALSE, xlab = '', ylab = '', xlim=xrange,
       ylim=yrange, cex.main=1.65, cex.axis=1.35, cex.lab=1.46,  main="", lwd=9, col='black')

  lines(xx, NAO.10yr.norm, type="l", lty = 2, cex=1.5, lwd=5.7, col='red3')
  xarr<- (0:7)*10+1950; numx<- length(xarr); yy1<- c(-1000,1000)
  for (x in 1:numx) { xx1<- c(xarr[x],xarr[x]); lines(xx1,yy1, type="l",lty = 2,lwd=1.65,col='grey') }
  yarr<- (-3:3)*1; numy<- length(yarr); xx1<- c(1900,2100)
  for (y in 1:numy) { yy1<- c(yarr[y],yarr[y]); lines(xx1,yy1, type="l",lty = 2,lwd=1.65,col='grey') }
  lines(c(1925,2045), c(0.0,0.0), type="l", lty = 1, cex=2.5, col='grey10')
  legend(1950, 2.5, c("Storm losses", "NAO fixed location"), lty=c(1,2), lwd=c(9,7), col=c('black', 'red3'), cex=1.75)

  ax.size = 1.7
  ## X Axis
  xticks <- seq(1950, 2020, by=10)
  axis(side=1, at = xticks, labels = xticks, las = 1, tck = 0.017, cex.axis = ax.size, padj=-0.26)
  ## Y Axis
  yticks <- (-6:6)*0.5;  axis(side=2, at = yticks, labels = NA, las = 1, tck = -0.0125, cex.axis = ax.size*1.15, hadj=0.65)
  yticks <- (-3:3)*1;  axis(side=2, at = yticks, labels = yticks, las = 1, tck = -0.0125, cex.axis = ax.size*1.15, hadj=0.65)
  mtext(side = 2, 'Standardised anomaly', cex = 3.75, padj = -2.9)
  yticks <- (-6:6)*0.5;  axis(side=4, at = yticks, labels = NA, las = 1, tck = 0.0125, cex.axis = ax.size*1.15, hadj=0.5)
box(); par(op); dev.off()



######################################################################################################################################################

###  Make map, probably not used in paper


init.map.plot <- function (fn.xlim = c(-10,30), fn.ylim = c(40,70), fn.xlab = "", fn.ylab = "", fn.cex.axis=1){
  plot(c(0,0), xlim = fn.xlim, ylim = fn.ylim, xaxs = 'i', yaxs = 'i',
       xlab = fn.xlab, ylab = fn.ylab, cex.axis= fn.cex.axis, add = FALSE) }

xlim1 <- -12; xlim2 <- 32; ylim1 <- 35; ylim2 <- 72

ctry.name.all<- c('Austria','Belgium','Switzerland','Germany','Denmark','France','Ireland','Luxembourg','Netherlands','Norway','Sweden','UK','Poland','Czech','Slovakia','Finland',
                  'Italy','Spain','Portugal','Andorra','Croatia','Bosnia','Montenegro','Albania','Russia','Slovenia','Serbia','Greece','North Macedonia','Romania',
				  'Latvia','Lithuania','Estonia','Ukraine','Belarus','Moldova','Cyprus','Turkey','Hungary','Algeria','Tunisia','Morocco','Libya','Egypt')

ctry.name.domain<- c('Austria','Belgium','Switzerland','Germany','Denmark','France','Ireland','Luxembourg','Netherlands','Norway','Sweden','UK','Poland','Czech','Slovakia','Finland')


op.fn<- paste0(opdir.plots, 'Figure_2_Map_backup.png')    #####  not used in submission due to restrictions in length of manuscript
png(file = op.fn, width = 1300., height = 900., type = 'cairo')
  op<- par( mai=c(0.45, 0.63, 0.3, 0.3), cex=1.9, mfrow=c(1,1))
  init.map.plot(fn.xlim = c(xlim1,xlim2), fn.ylim = c(ylim1,ylim2), fn.xlab = "", fn.ylab = "", fn.cex.axis=2.25)
  map('world', regions=ctry.name.all, add=TRUE, lwd=1.1, col='grey50')
  map('world', regions=ctry.name.domain, add=TRUE, lwd=1.65, fill=TRUE, col='grey83')
par(op); dev.off()





######################################################################################################################################################


### 2. compute annual loss index for each model and each country


# First, get the AAL (long-term mean) from the control run

sev.ev.ctl<- array(0.0, c(num.ctl.models, num.ctry));  freq.ev.ctl<- sev.ev.ctl;  aal.ctl<- sev.ev.ctl
ssi.rank.ctl<- array(0.0, c(num.ctl.models, num.ctry,6000));  aal.grid.ctl<- array(0.0, c(num.ctl.models, num.grid.pts))
ann.ctl.ctry<- array(NA, c(num.ctl.models, num.ctry, 600)); ann.ctl.ctry.freq<- array(0, c(num.ctl.models, num.ctry, 600))
ssi.ctl.thresh.keep<- array(0.0, num.ctl.models)
for (ff in 1:num.ctl.models) {

  ip.fn<- paste0(ipdir.ctl, model.id.ctl[ff], '_', ctl.runid.all[ff], '_', ctl.grid.all[ff], '_SSI_ev.csv')
  ctl.dat<- as.matrix(fread(ip.fn, header=F, sep=','));  row.names(ctl.dat)<- NULL ; colnames(ctl.dat)<- NULL
  num.days.per.yr<- dim(ctl.dat)[1]/num.yrs.ctl.model[ff]

  ssi.ev.ctl.tmp<- rowSums(ctl.dat)
  norm.fac<- 1e-6
  ssi.ev.ctl.tmp.2<- ssi.ev.ctl.tmp * norm.fac
  ssi.thresh<- 1e-3*(sum(ssi.ev.ctl.tmp.2)/num.yrs.ctl.model[ff])
  ssi.ctl.thresh.keep[ff]<- ssi.thresh
  ptr.ev.ctl<- which(ssi.ev.ctl.tmp.2 > ssi.thresh)
  storm.yr.num<- (ptr.ev.ctl %/% num.days.per.yr) + 1

  for (cc in 1:num.ctry) {
    ptr.ctry<- which(grid.ctry == all.ctry[cc]); if (cc == num.ctry) { ptr.ctry<- which(grid.ctry != 'XX') }

    ssi.ev.ctl<- rowSums(ctl.dat[,ptr.ctry])[ptr.ev.ctl] * norm.fac
    freq.ev.ctl[ff,cc]<- length(ptr.ev.ctl) / num.yrs.ctl.model[ff];  sev.ev.ctl[ff,cc]<- mean(ssi.ev.ctl[]);   aal.ctl[ff,cc]<- sum(ssi.ev.ctl)/num.yrs.ctl.model[ff]

    srt.ctl<- sort(ssi.ev.ctl, dec=T);  numpos<- min(6000, length(ptr.ev.ctl));  ssi.rank.ctl[ff, cc, 1:numpos]<- srt.ctl[1:numpos]

    for (g in 1:length(ptr.ctry)) {
      aal.grid.ctl[ff,ptr.ctry[g]]<- sum(ctl.dat[ptr.ev.ctl[], ptr.ctry[g]]) * norm.fac / num.yrs.ctl.model[ff]
    }
    for (yr in 1:num.yrs.ctl.model[ff]) {
      ptr.yr<- which(storm.yr.num == yr);  if (length(ptr.yr) > 0) { ann.ctl.ctry[ff, cc, yr]<- sum(ssi.ev.ctl[ptr.yr]) }
      thresh.temp<- ssi.ctl.thresh.keep[ff];  ptr.yr<- which(storm.yr.num == yr & ssi.ev.ctl > thresh.temp);  if (length(ptr.yr) > 0) { ann.ctl.ctry.freq[ff, cc, yr]<- length(ptr.yr) }
    }

  }
print(c(ff, substr(Sys.time(),12,19) ))
}


ann.ctl.ctry.filt<- ann.ctl.ctry; max.ctl.eu.filt<- array(0.0, num.ctl.models); min.ctl.eu.filt<- max.ctl.eu.filt; sd.ctl.eu<- max.ctl.eu.filt; mean.ctl.eu<- max.ctl.eu.filt; sd.ctl.eu.ann<- max.ctl.eu.filt
for (cc in 1:num.ctry) {
  for (ff in 1:num.ctl.models) {
    temparr<- 100.0 * (ann.ctl.ctry[ff, cc, 1:num.yrs.ctl.model[ff]]-aal.ctl[ff, cc]) / aal.ctl[ff, cc];  
    ann.ctl.ctry.filt[ff, cc, 1:num.yrs.ctl.model[ff]]<- pass.filt(temparr, filt.period, type="low", method="Butterworth", n=usr.filt.n)
    if (cc==num.ctry) {
      sd.ctl.eu.ann[ff]<- sd(temparr-ann.ctl.ctry.filt[ff, cc, 1:num.yrs.ctl.model[ff]])
      max.ctl.eu.filt[ff]<- max(ann.ctl.ctry.filt[ff, cc, 1:num.yrs.ctl.model[ff]])
      min.ctl.eu.filt[ff]<- min(ann.ctl.ctry.filt[ff, cc, 1:num.yrs.ctl.model[ff]])
      sd.ctl.eu[ff]<- sd(ann.ctl.ctry[ff, cc, 1:num.yrs.ctl.model[ff]])
      mean.ctl.eu[ff]<- mean(ann.ctl.ctry[ff, cc, 1:num.yrs.ctl.model[ff]])
    }
  }
}

###  print out standard errors of 20-yr mean storminess, relative to the mean storminess
for (ff in 1:num.ctl.models) {
  print(c( ff, sd.ctl.eu[ff], mean.ctl.eu[ff], (sd.ctl.eu[ff]/sqrt(20))/mean.ctl.eu[ff] ))
}

print( mean(sd.ctl.eu/sqrt(20)/mean.ctl.eu) )



max.anom.window<- c()
for (ff in 1:num.ctl.models) {
  num.windows<- trunc(num.yrs.ctl.model[ff]/72)
  ann.ctl.ctry.temp<- log(10*ann.ctl.ctry[ff, 17, 1:num.yrs.ctl.model[ff]])
  for (w in 1:num.windows) {
    indx<- (w-1)*72 + 1:72
    temparr<- ann.ctl.ctry.temp[indx]
    temparr2<- as.numeric(filter(temparr, rep(1/20, 20)))
    max.anom.window<- c(max.anom.window, max(temparr2,na.rm=T)-mean(temparr))
  }
}




############################################################

# Now, get the losses from the external forcing runs


ptr.CTL<- array(0, num.cmip6.models)
ptr.CTL[1:8] <- 1
ptr.CTL[9:17] <- 2
ptr.CTL[18:27] <- 3;  ptr.CTL[53:57] <- 3
ptr.CTL[28:37] <- 4
ptr.CTL[38:47] <- 5;  ptr.CTL[58:62] <- 5
ptr.CTL[48:52] <- 6

num.days.yr<- array(212,171); for (y in 0:42) { num.days.yr[(y*4+3)]<- 213 }; num.days.yr[51]<- 212
yr.hist<- c();  for (y in 1:171) { yr.hist<- c(yr.hist, rep(y, num.days.yr[y])) }

num.days.yr.360<- array(210,171)
yr.hist.360<- c();  for (y in 1:171) { yr.hist.360<- c(yr.hist.360, rep(y, num.days.yr.360[y])) }

num.days.yr.365<- array(213,171)
yr.hist.365<- c();  for (y in 1:171) { yr.hist.365<- c(yr.hist.365, rep(y, num.days.yr.365[y])) }

num.yrs.hist.model<- array(165, num.cmip6.models)
ptr.tmp<- c(38:47, 58:62)
num.yrs.damip.model<- array(171, num.cmip6.models); num.yrs.damip.model[ptr.tmp]<- 165


sev.ev.ghg<- array(0.0, c(num.cmip6.models, num.ctry)); freq.ev.ghg<- sev.ev.ghg; sev.ev.nat<- sev.ev.ghg; freq.ev.nat<- sev.ev.ghg; sev.ev.aero<- sev.ev.ghg; freq.ev.aero<- sev.ev.ghg
aal.ghg<- sev.ev.ghg; aal.nat<- sev.ev.ghg; aal.aero<- sev.ev.ghg
ssi.rank.ghg<- array(0.0, c(num.cmip6.models, num.ctry, 3000)); ssi.rank.nat<- ssi.rank.ghg; ssi.rank.aero<- ssi.rank.ghg
ann.ghg.ctry<- array(0.0, c(num.cmip6.models, num.ctry, 165)); ann.nat.ctry<- ann.ghg.ctry; ann.aero.ctry<- ann.ghg.ctry
aal.grid.ghg<- array(0.0, c(num.cmip6.models, num.grid.pts));  aal.grid.nat<- aal.grid.ghg; aal.grid.aero<- aal.grid.ghg
aal.grid.ghg.197099<- array(0.0, c(num.cmip6.models, num.grid.pts)); aal.grid.nat.197099<- aal.grid.ghg.197099; aal.grid.aero.197099<- aal.grid.ghg.197099
ann.ghg.ctry.freq<- array(0, c(num.cmip6.models, num.ctry, 165)); ann.nat.ctry.freq<- ann.ghg.ctry.freq; ann.aero.ctry.freq<- ann.ghg.ctry.freq

sev.ev.hist<- sev.ev.ghg; freq.ev.hist<- sev.ev.ghg; aal.hist<- sev.ev.ghg
ssi.rank.hist<- ssi.rank.ghg; ann.hist.ctry<- ann.ghg.ctry; aal.grid.hist<- aal.grid.ghg; aal.grid.hist.197099<- aal.grid.ghg.197099; 
ann.hist.ctry.freq<- array(0, c(num.cmip6.models, num.ctry, 165));  ann.hist.ctry<- array(0.0, c(num.cmip6.models, num.ctry, 165))

for (ff in 1:num.cmip6.models) {

  ip.fn<- paste0(ipdir.hist, model.id.all[ff], '_', hist.runid.all[ff], '_', hist.grid.all[ff], '_SSI_ev.csv')
  hist.dat<- as.matrix(fread(ip.fn, header=F, sep=','));  row.names(hist.dat)<- NULL ; colnames(hist.dat)<- NULL
  ip.fn<- paste0(ipdir.ghg, model.id.all[ff], '_', ghg.runid.all[ff], '_', ghg.grid.all[ff], '_SSI_ev.csv')
  ghg.dat<- as.matrix(fread(ip.fn, header=F, sep=','));  row.names(ghg.dat)<- NULL ; colnames(ghg.dat)<- NULL
  ip.fn<- paste0(ipdir.nat, model.id.all[ff], '_', nat.runid.all[ff], '_', nat.grid.all[ff], '_SSI_ev.csv')
  nat.dat<- as.matrix(fread(ip.fn, header=F, sep=','));  row.names(nat.dat)<- NULL ; colnames(nat.dat)<- NULL
  ip.fn<- paste0(ipdir.aero, model.id.all[ff], '_', aero.runid.all[ff], '_', aero.grid.all[ff], '_SSI_ev.csv')
  aero.dat<- as.matrix(fread(ip.fn, header=F, sep=','));  row.names(aero.dat)<- NULL ; colnames(aero.dat)<- NULL


  ssi.ev.ghg.tmp<- rowSums(ghg.dat); ssi.ev.nat.tmp<- rowSums(nat.dat); ssi.ev.aero.tmp<- rowSums(aero.dat); ssi.ev.hist.tmp<- rowSums(hist.dat)
  norm.fac<- 1e-6
  ssi.ev.ghg.tmp.2<- ssi.ev.ghg.tmp * norm.fac;  ssi.ev.nat.tmp.2<- ssi.ev.nat.tmp * norm.fac;  ssi.ev.aero.tmp.2<- ssi.ev.aero.tmp * norm.fac;  ssi.ev.hist.tmp.2<- ssi.ev.hist.tmp * norm.fac
  ssi.thresh<- ssi.ctl.thresh.keep[ptr.CTL[ff]]
# the following three variables are what's used in the next block of code, the loop over countries
  ptr.ev.ghg<- which(ssi.ev.ghg.tmp.2 > ssi.thresh); ptr.ev.nat<- which(ssi.ev.nat.tmp.2 > ssi.thresh); ptr.ev.aero<- which(ssi.ev.aero.tmp.2 > ssi.thresh); ptr.ev.hist<- which(ssi.ev.hist.tmp.2 > ssi.thresh)

  for (cc in 1:num.ctry) {
    ptr.ctry<- which(grid.ctry == all.ctry[cc]); if (cc == num.ctry) { ptr.ctry<- which(grid.ctry != 'XX') }

    ssi.ev.hist<- rowSums(hist.dat[,ptr.ctry])[ptr.ev.hist] * norm.fac
    freq.ev.hist[ff,cc]<- length(ptr.ev.hist) / num.yrs.hist.model[ff];  sev.ev.hist[ff,cc]<- mean(ssi.ev.hist[]);   aal.hist[ff,cc]<- sum(ssi.ev.hist)/num.yrs.hist.model[ff]

    ssi.ev.ghg<- rowSums(ghg.dat[,ptr.ctry])[ptr.ev.ghg] * norm.fac
    freq.ev.ghg[ff,cc]<- length(ptr.ev.ghg) / num.yrs.damip.model[ff];  sev.ev.ghg[ff,cc]<- mean(ssi.ev.ghg[]);   aal.ghg[ff,cc]<- sum(ssi.ev.ghg)/num.yrs.damip.model[ff]

    ssi.ev.nat<- rowSums(nat.dat[,ptr.ctry])[ptr.ev.nat] * norm.fac
    freq.ev.nat[ff,cc]<- length(ptr.ev.nat) / num.yrs.damip.model[ff];  sev.ev.nat[ff,cc]<- mean(ssi.ev.nat[]);   aal.nat[ff,cc]<- sum(ssi.ev.nat)/num.yrs.damip.model[ff]

    ssi.ev.aero<- rowSums(aero.dat[,ptr.ctry])[ptr.ev.aero] * norm.fac
    freq.ev.aero[ff,cc]<- length(ptr.ev.aero) / num.yrs.damip.model[ff];  sev.ev.aero[ff,cc]<- mean(ssi.ev.aero[]);   aal.aero[ff,cc]<- sum(ssi.ev.aero)/num.yrs.damip.model[ff]

    srt.hist<- sort(ssi.ev.hist, dec=T);  numpos<- min(3000, length(ptr.ev.hist));  ssi.rank.hist[ff, cc, 1:numpos]<- srt.hist[1:numpos]
    srt.ghg<- sort(ssi.ev.ghg, dec=T);  numpos<- min(3000, length(ptr.ev.ghg));  ssi.rank.ghg[ff, cc, 1:numpos]<- srt.ghg[1:numpos]
    srt.nat<- sort(ssi.ev.nat, dec=T);  numpos<- min(3000, length(ptr.ev.nat));  ssi.rank.nat[ff, cc, 1:numpos]<- srt.nat[1:numpos]
    srt.aero<- sort(ssi.ev.aero, dec=T);  numpos<- min(3000, length(ptr.ev.aero));  ssi.rank.aero[ff, cc, 1:numpos]<- srt.aero[1:numpos]

    for (g in 1:length(ptr.ctry)) {
      aal.grid.hist[ff,ptr.ctry[g]]<- sum(hist.dat[ptr.ev.hist[], ptr.ctry[g]]) * norm.fac / num.yrs.hist.model[ff]
      aal.grid.ghg[ff,ptr.ctry[g]]<- sum(ghg.dat[ptr.ev.ghg[], ptr.ctry[g]]) * norm.fac / num.yrs.damip.model[ff]
      aal.grid.nat[ff,ptr.ctry[g]]<- sum(nat.dat[ptr.ev.nat[], ptr.ctry[g]]) * norm.fac / num.yrs.damip.model[ff]
      aal.grid.aero[ff,ptr.ctry[g]]<- sum(aero.dat[ptr.ev.aero[], ptr.ctry[g]]) * norm.fac / num.yrs.damip.model[ff]
    }

    if (model.id.all[ff] == "HadGEM3-GC31-MM" | model.id.all[ff] == "HadGEM3-GC31-LL" | model.id.all[ff] == "KACE-1-0-G" | model.id.all[ff] == "UKESM1-0-LL") {
      yrs.arr<- yr.hist.360
    } else if (model.id.all[ff] == "CanESM5" | model.id.all[ff] == "CMCC-CM2-SR5" | model.id.all[ff] == "GFDL-ESM4") {
      yrs.arr<- yr.hist.365
    } else { yrs.arr<- yr.hist }
    for (yr in 1:165) {
      ptr.yr1<- which(yrs.arr[ptr.ev.hist] == yr);  if (length(ptr.yr1) > 0) { ann.hist.ctry[ff,cc,yr]<- sum(ssi.ev.hist[ptr.yr1]) }
      ptr.yr2<- which(yrs.arr[ptr.ev.ghg] == yr);  if (length(ptr.yr2) > 0) { ann.ghg.ctry[ff,cc,yr]<- sum(ssi.ev.ghg[ptr.yr2]) }
      ptr.yr3<- which(yrs.arr[ptr.ev.nat] == yr);  if (length(ptr.yr3) > 0) { ann.nat.ctry[ff,cc,yr]<- sum(ssi.ev.nat[ptr.yr3]) }
      ptr.yr4<- which(yrs.arr[ptr.ev.aero] == yr);  if (length(ptr.yr4) > 0) { ann.aero.ctry[ff,cc,yr]<- sum(ssi.ev.aero[ptr.yr4]) }
      thresh.temp<- ssi.ctl.thresh.keep[ptr.CTL[ff]];  ptr.yr1<- which(yrs.arr[ptr.ev.hist] == yr & ssi.ev.hist > thresh.temp);  if (length(ptr.yr1) > 0) { ann.hist.ctry.freq[ff,cc,yr]<- length(ptr.yr1) }
      thresh.temp<- ssi.ctl.thresh.keep[ptr.CTL[ff]];  ptr.yr2<- which(yrs.arr[ptr.ev.ghg] == yr & ssi.ev.ghg > thresh.temp);  if (length(ptr.yr2) > 0) { ann.ghg.ctry.freq[ff,cc,yr]<- length(ptr.yr2) }
      thresh.temp<- ssi.ctl.thresh.keep[ptr.CTL[ff]];  ptr.yr3<- which(yrs.arr[ptr.ev.nat] == yr & ssi.ev.nat > thresh.temp);  if (length(ptr.yr3) > 0) { ann.nat.ctry.freq[ff,cc,yr]<- length(ptr.yr3) }
      thresh.temp<- ssi.ctl.thresh.keep[ptr.CTL[ff]];  ptr.yr4<- which(yrs.arr[ptr.ev.aero] == yr & ssi.ev.aero > thresh.temp);  if (length(ptr.yr4) > 0) { ann.aero.ctry.freq[ff,cc,yr]<- length(ptr.yr4) }
    }
	
    ptr.yr1<- which(yrs.arr[ptr.ev.hist] >= 121 & yrs.arr[ptr.ev.hist] <= 150)
    ptr.yr2<- which(yrs.arr[ptr.ev.ghg] >= 121 & yrs.arr[ptr.ev.ghg] <= 150)
    ptr.yr3<- which(yrs.arr[ptr.ev.nat] >= 121 & yrs.arr[ptr.ev.nat] <= 150)
    ptr.yr4<- which(yrs.arr[ptr.ev.aero] >= 121 & yrs.arr[ptr.ev.aero] <= 150)
    for (g in 1:length(ptr.ctry)) {
      aal.grid.hist.197099[ff,ptr.ctry[g]]<- sum(hist.dat[ptr.ev.hist[ptr.yr1], ptr.ctry[g]]) * norm.fac / 30
      aal.grid.ghg.197099[ff,ptr.ctry[g]]<- sum(ghg.dat[ptr.ev.ghg[ptr.yr2], ptr.ctry[g]]) * norm.fac / 30
      aal.grid.nat.197099[ff,ptr.ctry[g]]<- sum(nat.dat[ptr.ev.nat[ptr.yr3], ptr.ctry[g]]) * norm.fac / 30
      aal.grid.aero.197099[ff,ptr.ctry[g]]<- sum(aero.dat[ptr.ev.aero[ptr.yr4], ptr.ctry[g]]) * norm.fac / 30
    }
	
  }
print(c(ff, substr(Sys.time(),12,19) ))
}


ann.ghg.eu.freq<- array(NA, c(num.cmip6.models, 165));  ann.nat.eu.freq<- ann.ghg.eu.freq;  ann.aero.eu.freq<- ann.ghg.eu.freq
for (ff in 1:num.cmip6.models) {
  ann.ghg.eu.freq[ff, ]<- 100.0 * (ann.ghg.ctry.freq[ff, 17, ]-freq.ev.ctl[ptr.CTL[ff], 17]) / freq.ev.ctl[ptr.CTL[ff], 17]
  ann.nat.eu.freq[ff, ]<- 100.0 * (ann.nat.ctry.freq[ff, 17, ]-freq.ev.ctl[ptr.CTL[ff], 17]) / freq.ev.ctl[ptr.CTL[ff], 17]
  ann.aero.eu.freq[ff, ]<- 100.0 * (ann.aero.ctry.freq[ff, 17, ]-freq.ev.ctl[ptr.CTL[ff], 17]) / freq.ev.ctl[ptr.CTL[ff], 17]
}

ann.ghg.eu.freq.mn<- colMeans(ann.ghg.eu.freq, na.rm=T)
ann.nat.eu.freq.mn<- colMeans(ann.nat.eu.freq, na.rm=T)
ann.aero.eu.freq.mn<- colMeans(ann.aero.eu.freq, na.rm=T)



ann.ghg.eu.aal<- array(NA, c(num.cmip6.models, 165));  ann.nat.eu.aal<- ann.ghg.eu.aal;  ann.aero.eu.aal<- ann.ghg.eu.aal
for (ff in 1:num.cmip6.models) {
  ann.ghg.eu.aal[ff, ]<- 100.0 * (ann.ghg.ctry[ff, 17, ]-aal.ctl[ptr.CTL[ff], 17]) / aal.ctl[ptr.CTL[ff], 17]
  ann.nat.eu.aal[ff, ]<- 100.0 * (ann.nat.ctry[ff, 17, ]-aal.ctl[ptr.CTL[ff], 17]) / aal.ctl[ptr.CTL[ff], 17]
  ann.aero.eu.aal[ff, ]<- 100.0 * (ann.aero.ctry[ff, 17, ]-aal.ctl[ptr.CTL[ff], 17]) / aal.ctl[ptr.CTL[ff], 17]
}

ann.ghg.eu.aal.mn<- colMeans(ann.ghg.eu.aal, na.rm=T)
ann.nat.eu.aal.mn<- colMeans(ann.nat.eu.aal, na.rm=T)
ann.aero.eu.aal.mn<- colMeans(ann.aero.eu.aal, na.rm=T)



ann.ghg.ctry.filt<- array(NA, c(num.cmip6.models, num.ctry, 165));  ann.nat.ctry.filt<- ann.nat.ctry;  ann.aero.ctry.filt<- ann.aero.ctry
for (cc in 1:num.ctry) {
  for (ff in 1:num.cmip6.models) {
    temparr<- 100.0 * (ann.ghg.ctry[ff, cc, ]-aal.ctl[ptr.CTL[ff], cc]) / aal.ctl[ptr.CTL[ff], cc]; ann.ghg.ctry.filt[ff, cc, ]<- pass.filt(temparr, filt.period, type="low", method="Butterworth", n=usr.filt.n)
    temparr<- 100.0 * (ann.nat.ctry[ff, cc, ]-aal.ctl[ptr.CTL[ff], cc]) / aal.ctl[ptr.CTL[ff], cc]; ann.nat.ctry.filt[ff, cc, ]<- pass.filt(temparr, filt.period, type="low", method="Butterworth", n=usr.filt.n)
    temparr<- 100.0 * (ann.aero.ctry[ff, cc, ]-aal.ctl[ptr.CTL[ff], cc]) / aal.ctl[ptr.CTL[ff], cc]; ann.aero.ctry.filt[ff, cc, ]<- pass.filt(temparr, filt.period, type="low", method="Butterworth", n=usr.filt.n)
  }
}



usr.ind<- c(1:num.cmip6.models)
aal.ghg.ensmn<- colMeans(ann.ghg.ctry.filt[usr.ind, , ], na.rm=T)
aal.nat.ensmn<- colMeans(ann.nat.ctry.filt[usr.ind, , ], na.rm=T)
aal.aero.ensmn<- colMeans(ann.aero.ctry.filt[usr.ind, , ], na.rm=T)

aal.ghg.ensmn.eu<- aal.ghg.ensmn[17, ]
aal.nat.ensmn.eu<- aal.nat.ensmn[17, ]
aal.aero.ensmn.eu<- aal.aero.ensmn[17, ]


ann.aero.ctry.testing<- array(NA, c(num.cmip6.models, 3))
for (ff in 1:num.cmip6.models) {
  temparr<- 100.0 * (ann.aero.ctry[ff, 17, ]-aal.ctl[ptr.CTL[ff], 17]) / aal.ctl[ptr.CTL[ff], 17]
  ann.aero.ctry.testing[ff,1]<- mean(temparr[111:130]); ann.aero.ctry.testing[ff,2]<- mean(temparr[131:150]); ann.aero.ctry.testing[ff,3]<- mean(temparr[151:165])
}
colMeans(ann.aero.ctry.testing)




uniq.models<- unique(model.id.all); num.uniq.models<- length(uniq.models)

ann.ghg.ctry.filt.model<- array(NA, c(num.uniq.models, num.ctry, 165)); ann.nat.ctry.filt.model<- ann.ghg.ctry.filt.model; ann.aero.ctry.filt.model<- ann.ghg.ctry.filt.model
for (mm in 1:num.uniq.models) {
  mod.indx<- which(model.id.all == uniq.models[mm])
  temparr<- colMeans(ann.ghg.ctry.filt[mod.indx, , ], na.rm=T);  for (cc in 1:num.ctry) { for (yy in 1:165) { ann.ghg.ctry.filt.model[mm, cc, yy]<- temparr[cc,yy] } }
  temparr<- colMeans(ann.nat.ctry.filt[mod.indx, , ], na.rm=T);  for (cc in 1:num.ctry) { for (yy in 1:165) { ann.nat.ctry.filt.model[mm, cc, yy]<- temparr[cc,yy] } }
  temparr<- colMeans(ann.aero.ctry.filt[mod.indx, , ], na.rm=T);  for (cc in 1:num.ctry) { for (yy in 1:165) { ann.aero.ctry.filt.model[mm, cc, yy]<- temparr[cc,yy] } }
}

aal.ghg.model.eu<- array(NA, c(num.uniq.models, 165));  aal.nat.model.eu<- aal.ghg.model.eu;  aal.aero.model.eu<- aal.ghg.model.eu
for (mm in 1:num.uniq.models) {
  aal.ghg.model.eu[mm, 1:165]<- ann.ghg.ctry.filt.model[mm, 17, 1:165]
  aal.nat.model.eu[mm, 1:165]<- ann.nat.ctry.filt.model[mm, 17, 1:165]
  aal.aero.model.eu[mm, 1:165]<- ann.aero.ctry.filt.model[mm, 17, 1:165]
}






#######################################################################################################################################

### Process, and make Figure 4 on internal variability in Control runs


###  compute the size of the peak 20-year mean loss, versus average loss in 75-year chunks from control runs
max.anom.window<- c()
for (ff in 1:num.ctl.models) {
  num.windows<- trunc(num.yrs.ctl.model[ff]/75)
  ann.ctl.ctry.temp<- ann.ctl.ctry[ff, 17, 1:num.yrs.ctl.model[ff]]
  for (w in 1:num.windows) {
    indx<- (w-1)*75 + 1:75
    temparr<- ann.ctl.ctry.temp[indx]
    temparr2<- as.numeric(filter(temparr, rep(1/20, 20)))
    max.anom.window<- c(max.anom.window, max(temparr2,na.rm=T)/mean(temparr))
  }
}

###  compute the same as above, from Aero history runs
aero.max.anom.window<- array(NA, num.cmip6.models)
for (ff in 1:num.cmip6.models) {
  temparr<- ann.aero.ctry[ff, 17, 91:165]
  temparr2 <- as.numeric(filter(temparr, rep(1/20, 20)))
  aero.max.anom.window[ff] <- max(temparr2, na.rm=T) / mean(temparr)
}


### calculate the histogram data, for the plot
usr.breks<- ((1:9)*0.1)+1.0; num.bins<- length(usr.breks)-1
num.occur.bins.ctl<- array(0, num.bins);  num.occur.bins.aero<- num.occur.bins.ctl
for (b in 1:num.bins) {
  num.occur.bins.ctl[b]<- length(which(max.anom.window > usr.breks[b] & max.anom.window <= usr.breks[(b+1)]))
  num.occur.bins.aero[b]<- length(which(aero.max.anom.window > usr.breks[b] & aero.max.anom.window <= usr.breks[(b+1)]))
}
freq.occur.bins.ctl<- num.occur.bins.ctl / sum(num.occur.bins.ctl)
freq.occur.bins.aero<- num.occur.bins.aero / sum(num.occur.bins.aero)



###  compute the same from observed

ev.loss.fn<- 'D:/***/**.csv'
ip.dat<- read.table(ev.loss.fn, header=T, sep=',')
px<- which(as.integer(ip.dat$Year) >= 1940 & as.integer(ip.dat$Year) <= 2014)
ip.date<- as.integer(ip.dat$Year[px])
euws.ann.loss<- as.numeric(ip.dat$EU[px])*1e-6
obs.max.anom<- mean(euws.ann.loss[41:60]) / mean(euws.ann.loss)




fn<- paste0(opdir.plots, "Figure_4.png")
png(file = fn, width = 1300., height = 700., res=216, pointsize=10 )
op <- par( mai = c(0.55, 0.67, 0.1, 0.1), cex=1.0, mfrow=c(1,1))
ymax.val<- 0.4; ymin.val<- 0.0; yrange<- c(ymin.val,ymax.val+0.02)

barplot(freq.occur.bins.aero, space=0.0, col='lightgreen', ylim=yrange, main='', las=1, names.arg=c(NA,NA,NA,NA,NA,NA,NA,NA), ylab='', 
                                                                     cex.axis=1.1, cex.main=1.3, cex.names=1.0, cex.lab=1.3, axes=F) 
xx<- c(-50,50); yarr<- seq(0.1,ymax.val, by=0.1); numy<- length(yarr)
for (y in 1:numy) {yy<- c(yarr[y],yarr[y]); lines(xx,yy,type="l",lty = 2,lwd=1,cex=1.0,col='grey30')}

barplot(freq.occur.bins.aero, space=0.0, col='lightgreen', ylim=yrange, main='', las=1, names.arg=c(NA,NA,NA,NA,NA,NA,NA,NA), ylab='', 
                                                                     cex.axis=1.1, cex.main=1.3, cex.names=1.0, cex.lab=1.3, axes=F, add=TRUE) 
obs.xval<- (obs.max.anom-1.1)/0.8 * 8
xx<- c(obs.xval, obs.xval);  yy<-  c(0, 0.35);  lines(xx, yy, type="l", lty = 1, cex=1.5, lwd=4.0, col='black')
  xticks <- seq(0.0, 8.0, by=1);  xtick.lbl<- seq(1.1, 1.9, by=0.1)
  axis(side=1, at = xticks, labels = xtick.lbl, las = 1, tck = -0.0135, cex.axis = 1.15, padj=-0.91)
  yticks <- seq(ymin.val,ymax.val, by=0.1)
  axis(side=2, at = yticks, labels = yticks, las = 1, tck = -0.015, cex.axis = 1.15, hadj=0.475)
  mtext(side = 2, 'Probability mass', cex = 1.3, padj = -3.25)
  mtext(side = 1, 'Maximum 20-year  loss  / average loss', cex = 1.2, padj = 2.45)
#text(0.35, 0.36, '(a)', cex=1.95, col='black')
box(); par(op); dev.off()





#######################################################################################################################################

###  Make Figure 2 of article



### Make standard errors for the multi-model mean, and each model mean - used in Figure 2


ann.ghg.eu.all<- array(0.0, c(num.cmip6.models, 165));  ann.nat.eu.all<- ann.ghg.eu.all;  ann.aero.eu.all<- ann.ghg.eu.all
for (ff in 1:num.cmip6.models) {
  ann.ghg.eu.all[ff,]<- ann.ghg.ctry.filt[ff, 17, ]*aal.ctl[ptr.CTL[ff], 17]*0.01  ### 
  ann.nat.eu.all[ff,]<- ann.nat.ctry.filt[ff, 17, ]*aal.ctl[ptr.CTL[ff], 17]*0.01  ### 
  ann.aero.eu.all[ff,]<- ann.aero.ctry.filt[ff, 17, ]*aal.ctl[ptr.CTL[ff], 17]*0.01  ### 
}

stdev.ghg.member<- array(0.0, 165);  stdev.nat.member<- stdev.ghg.member;  stdev.aero.member<- stdev.ghg.member
for (yr in 1:165) {
  stdev.ghg.member[yr]<- 100 * sd(ann.ghg.eu.all[ , yr]) / mean(aal.ctl[, 17])
  stdev.nat.member[yr]<- 100 * sd(ann.nat.eu.all[ , yr]) / mean(aal.ctl[, 17])
  stdev.aero.member[yr]<- 100 * sd(ann.aero.eu.all[ , yr]) / mean(aal.ctl[, 17])
}

stderr.ghg.ensmn<- stdev.ghg.member / sqrt(num.cmip6.models)
stderr.nat.ensmn<- stdev.nat.member / sqrt(num.cmip6.models)
stderr.aero.ensmn<- stdev.aero.member / sqrt(num.cmip6.models)


uniq.models<- unique(model.id.all); num.uniq.models<- length(uniq.models)

stderr.ghg.model.mn<- array(0.0, c(num.uniq.models, 165)); stderr.nat.model.mn<- stderr.ghg.model.mn; stderr.aero.model.mn<- stderr.ghg.model.mn
for (mm in 1:num.uniq.models) {
  mod.indx<- which(model.id.all == uniq.models[mm])
  for (yr in 1:165) {
    temparr<- 100 * sd(ann.ghg.eu.all[mod.indx, yr]) / mean(aal.ctl[mm, 17])   ###  std dev in %
    stderr.ghg.model.mn[mm, yr]<- temparr / sqrt(length(mod.indx))
    temparr<- 100 * sd(ann.nat.eu.all[mod.indx, yr]) / mean(aal.ctl[mm, 17])   ###  std dev in %
    stderr.nat.model.mn[mm, yr]<- temparr / sqrt(length(mod.indx))
    temparr<- 100 * sd(ann.aero.eu.all[mod.indx, yr]) / mean(aal.ctl[mm, 17])   ###  std dev in %
    stderr.aero.model.mn[mm, yr]<- temparr / sqrt(length(mod.indx))
 }
}

ann.aero.ctry.diff<- array(NA, c(num.cmip6.models, 165))
for (ff in 1:num.cmip6.models) {
  ann.aero.ctry.diff[ff, ]<- 100.0 * (ann.aero.ctry[ff, 17, ]-aal.ctl[ptr.CTL[ff], 17]) / aal.ctl[ptr.CTL[ff], 17]
}

ann.aero.ctry.model.diff<- array(NA, c(num.uniq.models, 165))
for (mm in 1:num.uniq.models) {
  mod.indx<- which(model.id.all == uniq.models[mm])
  temparr<- colMeans(ann.aero.ctry.diff[mod.indx, ], na.rm=T);  for (yy in 1:165) { ann.aero.ctry.model.diff[mm, yy]<- temparr[yy] }
}

model.pval<- array(0.0, num.uniq.models)
for (mm in 1:num.uniq.models) {   zz<- t.test(ann.aero.ctry.model.diff[mm, 121:160]);  model.pval[mm]<- 1e-4*round(1e4*zz$p.val) }



### plot out Figure 2
 
fn<- paste0(opdir.plots, "Figure_2.png")
png(file = fn, width = 2000., height = 2000., res=216, pointsize=10, type = 'cairo')
#op <- par( mai = c(0.4, 0.85, 0.3, 0.2), cex=1.0, mfrow=c(1,1)); yrs.plt<- 1850:2014; yrs.plt.obs<- 1940:2022
op <- par( mai = c(0.53, 0.77, 0.1, 0.1), cex=1.0, mfrow=c(2,1)); yrs.plt<- 1850:2014; yrs.plt.obs<- 1940:2022

  yy.max<- 65.0; yy.min<- -20;  del.yy<- 20.0
  yrange=c(yy.min-5, yy.max); xrange=c(1851,2014)
  plot(yrs.plt, aal.aero.ensmn.eu, type="l", lty = 1, cex=2.5, axes = FALSE, xlab = '', ylab = '', xlim=xrange,
       ylim=yrange, cex.main=1.65, cex.axis=1.35, cex.lab=1.46, main="", lwd=4.0, col='blue3')
  zz1<- aal.aero.ensmn.eu + 2*stderr.aero.ensmn; lines(yrs.plt, zz1, type="l", lty = 2, cex=1.5, lwd=1.15, col='dodgerblue2')
  zz2<- aal.aero.ensmn.eu - 2*stderr.aero.ensmn; lines(yrs.plt, zz2, type="l", lty = 2, cex=1.5, lwd=1.15, col='dodgerblue2')

  lines(yrs.plt, aal.aero.ensmn.eu, type="l", lty = 1, cex=1.5, lwd=4.0, col='blue3')

  lines( c(1800,2100), c(0,0), type="l", lty = 1, lwd=1.5, col='black')
  xarr<- (0:9)*20+1850; numx<- length(xarr); yy1<- c(-200,200)
  for (x in 1:numx) { xx1<- c(xarr[x],xarr[x]); lines(xx1,yy1, type="l",lty = 3,lwd=1.03,col='grey') }
  yarr<- (-8:8)*20; numy<- length(yarr); xx1<- c(1800,2100)
  for (y in 1:numy) { yy1<- c(yarr[y],yarr[y]); lines(xx1,yy1, type="l",lty = 3,lwd=1.03,col='grey') }

  text(1855, 54.5, '(a)', cex=2.05, col='black')
  ax.size = 1.65
  ## X Axis
  xticks <- seq(1850, 2010, by=20)
  axis(side=1, at = xticks, labels = xticks, las = 1, tck = 0.015, cex.axis = ax.size, padj=-0.55)
  ## Y Axis
  yticks <- (-8:8)*20
  axis(side=2, at = yticks, labels = yticks, las = 1, tck = -0.012, cex.axis = ax.size, hadj=0.7)
  mtext(side = 2, 'change in annual loss (%)', cex = 1.7, padj = -3.0)
box()


  rgb.palette2 <- colorRampPalette(c( "blue4","blue1","lightblue4"), space ="rgb")
  plotclr <- c('black', rgb.palette2(4),'darkolivegreen3')

  yy.max<- 130.0; yy.min<- -30;  del.yy<- 25.0
  yrange=c(yy.min, yy.max); xrange=c(1851,2014)
  plot(yrs.plt, aal.aero.ensmn.eu, type="l", lty = 0, cex=2.5, axes = FALSE, xlab = '', ylab = '', xlim=xrange,
       ylim=yrange, cex.main=1.65, cex.axis=1.35, cex.lab=1.46, main="", lwd=0, col='blue3')
  for (mm in 1:num.uniq.models) {  temparr<- aal.aero.model.eu[mm, ]; lines(yrs.plt, temparr, type="l", lty = mm, cex=1.3, lwd=2.3, col=plotclr[mm])  }
  lines( c(1800,2100), c(0,0), type="l", lty = 1, lwd=1.5, col='black')
  xarr<- (0:9)*20+1850; numx<- length(xarr); yy1<- c(-200,200)
  for (x in 1:numx) { xx1<- c(xarr[x],xarr[x]); lines(xx1,yy1, type="l",lty = 3,lwd=1.03,col='grey') }
  yarr<- (-8:8)*del.yy; numy<- length(yarr); xx1<- c(1800,2100)
  for (y in 1:numy) { yy1<- c(yarr[y],yarr[y]); lines(xx1,yy1, type="l",lty = 3,lwd=1.03,col='grey') }
  temparr<- sprintf("%1.4f", model.pval)
  legend(1867, yy.max+2, c(paste0("CMCC-CM2-SR5  (",temparr[1],")"), paste0("CanESM5  (",temparr[2],")"), paste0("HadGEM3-GC31-LL  (",temparr[3],")"), paste0("MIROC6  (",temparr[4],")"), 
                           paste0("MPI-ESM1-2-LR  (",temparr[5],")"), paste0("MRI-ESM2-0  (",temparr[6],")")),  lty=c(1:6), lwd=c(3,3,3,3,3,3), pch=NA, col=c(plotclr), cex=1.38)
  text(1855, 116, '(b)', cex=2.05, col='black')
  ax.size = 1.65
  ## X Axis
  xticks <- seq(1850, 2010, by=20)
  axis(side=1, at = xticks, labels = xticks, las = 1, tck = 0.015, cex.axis = ax.size, padj=-0.55)
  ## Y Axis
  yticks <- (-8:8)*del.yy
  axis(side=2, at = yticks, labels = yticks, las = 1, tck = -0.012, cex.axis = ax.size, hadj=0.7)
  mtext(side = 2, 'change in annual loss (%)', cex = 1.7, padj = -3.12)
box(); par(op); dev.off()





################################################################################################################################################################

###  analyse surface pressure and thicknesses


ipdir.ip<- "D:/***/"
ipdir.gz.ctl<- "E:/***/"

ip.fn<- paste0(ipdir.ip, 'CMIP6_model_names_DAMIP_ps_gz_CTL.csv')
ip.data<- read.table(ip.fn, header=T, sep=',')
model.id.ctl<- as.character(ip.data$model); ctl.runid.all<- as.character(ip.data$CTL_run); ctl.grid.all<- as.character(ip.data$CTL_grid)
num.ctl.models<- length(model.id.ctl)

num.yrs.ctl.model<- c(500, 600, 500, 500, 600, 200)

ip.fn<- paste0(ipdir.ip, 'universal_grid_1p0.csv')
ip.data<- read.table(ip.fn, header=T, sep=',')
grid.lat.psgz<- ip.data$LatC;  grid.lon.psgz<- ip.data$LonC
num.grid.pts.psgz<- length(grid.lat.psgz)
grid.wt.psgz<- cos(grid.lat.psgz * pi / 180)


ctl.gz.mn<- array(0.0, c(num.ctl.models, num.grid.pts.psgz));  ctl.ps.mn<- ctl.gz.mn
for (ff in 1:num.ctl.models) {

  ip.fn<- paste0(ipdir.gz.ctl, model.id.ctl[ff], '_z500_ann.csv')
  ctl.dat<- fread(ip.fn, header=F, sep=',')
  ctl.gz.mn[ff, 1:num.grid.pts.psgz] <- colMeans(ctl.dat)

  ip.fn2<- paste0(ipdir.gz.ctl, model.id.ctl[ff], '_ps_ann.csv')
  ctl.dat2<- fread(ip.fn2, header=F, sep=',')
  ctl.ps.mn[ff, 1:num.grid.pts.psgz] <- colMeans(ctl.dat2)

print(c(ff, substr(Sys.time(),12,19) ))
}




ipdir.int<- "E:/***/"
ipdir.tabi.aero<- paste0(ipdir.int, "HIST_Aer/global_grid/")
ipdir.tabi.ghg<- paste0(ipdir.int, "HIST_GHG/global_grid/")
ipdir.tabi.nat<- paste0(ipdir.int, "HIST_Nat/global_grid/")

ip.fn<- paste0(ipdir.ip, 'CMIP6_model_names_DAMIP_ps_gz.csv')
ip.data<- read.table(ip.fn, header=T, sep=',')
model.id.test<- as.character(ip.data$model); aero.runid.test<- as.character(ip.data$Hist_Aer_run); aero.grid.test<- as.character(ip.data$Hist_Aer_grid)
ghg.runid.test<- as.character(ip.data$Hist_GHG_run); ghg.grid.test<- as.character(ip.data$Hist_GHG_grid)
nat.runid.test<- as.character(ip.data$Hist_Nat_run); nat.grid.test<- as.character(ip.data$Hist_Nat_grid)

num.test.models<- length(model.id.test)

sim.yrs<- 1850:2014;  num.winters<- length(sim.yrs)

test.yrs<- 121:160  ##  1970-2009
test.yrs.1940<- 91:165  ##  1940-2014
test.yrs.nat<- 143:147  ##  1992-1996


aero.ps.all<-array(0.0, c(num.test.models, num.grid.pts.psgz));  aero.gz.all<- aero.ps.all;  ghg.ps.all<- aero.ps.all;  ghg.gz.all<- aero.ps.all;  nat.ps.all<- aero.ps.all;  nat.gz.all<- aero.ps.all
aero.gz.all.1940<- aero.ps.all; ghg.gz.all.1940<- aero.ps.all; nat.gz.all.1940<- aero.ps.all
nat.gz.all.ann<- array(0.0, c(num.test.models, 30, num.grid.pts.psgz))
for (ff in 1:num.test.models) {
  ip.fn2<- paste0(ipdir.tabi.aero, model.id.test[ff], '_', aero.runid.test[ff], '_', aero.grid.test[ff], '_z500_ann.csv')
  ip.dat2<- fread(ip.fn2, header=F, sep=',')
  aero.gz.all[ff, 1:num.grid.pts.psgz] <- colMeans(ip.dat2[test.yrs,])
  aero.gz.all.1940[ff, 1:num.grid.pts.psgz] <- colMeans(ip.dat2[test.yrs.1940,])

  ip.fn2<- paste0(ipdir.tabi.ghg, model.id.test[ff], '_', ghg.runid.test[ff], '_', ghg.grid.test[ff], '_z500_ann.csv')
  ip.dat2<- fread(ip.fn2, header=F, sep=',')
  ghg.gz.all[ff, 1:num.grid.pts.psgz] <- colMeans(ip.dat2[test.yrs,])
  ghg.gz.all.1940[ff, 1:num.grid.pts.psgz] <- colMeans(ip.dat2[test.yrs.1940,])

  ip.fn2<- paste0(ipdir.tabi.nat, model.id.test[ff], '_', nat.runid.test[ff], '_', nat.grid.test[ff], '_z500_ann.csv')
  ip.dat2<- fread(ip.fn2, header=F, sep=',')
  nat.gz.all[ff, 1:num.grid.pts.psgz] <- colMeans(ip.dat2[test.yrs.nat,])
  nat.gz.all.1940[ff, 1:num.grid.pts.psgz] <- colMeans(ip.dat2[test.yrs.1940,])
  for (y in 136:165) { y2<- y-135; nat.gz.all.ann[ff, y2, 1:num.grid.pts.psgz] <- as.numeric(ip.dat2[y, 1:num.grid.pts.psgz]) }

  ip.fn<- paste0(ipdir.tabi.aero, model.id.test[ff], '_', aero.runid.test[ff], '_', aero.grid.test[ff], '_ps_ann.csv')
  ip.dat<- fread(ip.fn, header=F, sep=',')
  aero.ps.all[ff, 1:num.grid.pts.psgz] <- colMeans(ip.dat[test.yrs,])

  ip.fn<- paste0(ipdir.tabi.ghg, model.id.test[ff], '_', ghg.runid.test[ff], '_', ghg.grid.test[ff], '_ps_ann.csv')
  ip.dat<- fread(ip.fn, header=F, sep=',')
  ghg.ps.all[ff, 1:num.grid.pts.psgz] <- colMeans(ip.dat[test.yrs,])

  ip.fn<- paste0(ipdir.tabi.nat, model.id.test[ff], '_', nat.runid.test[ff], '_', nat.grid.test[ff], '_ps_ann.csv')
  ip.dat<- fread(ip.fn, header=F, sep=',')
  nat.ps.all[ff, 1:num.grid.pts.psgz] <- colMeans(ip.dat[test.yrs.nat,])

print(c(ff, substr(Sys.time(),12,19) ))
}



ptr.CTL<- array(0, num.test.models)
ptr.CTL[1:8] <- 1
ptr.CTL[9:17] <- 2
ptr.CTL[18:27] <- 3;  ptr.CTL[53:57] <- 3
ptr.CTL[28:37] <- 4
ptr.CTL[38:47] <- 5;  ptr.CTL[58:62] <- 5
ptr.CTL[48:52] <- 6


aero.ps.anom<- array(0.0, c(num.test.models, num.grid.pts.psgz)); aero.gz.anom<- aero.ps.anom; ghg.ps.anom<- aero.ps.anom; ghg.gz.anom<- aero.ps.anom; nat.ps.anom<- aero.ps.anom; nat.gz.anom<- aero.ps.anom
aero.gz.anom.1940<- aero.ps.anom;  ghg.gz.anom.1940<- aero.ps.anom;  nat.gz.anom.1940<- aero.ps.anom
nat.gz.anom.ann<- array(0.0, c(num.test.models, 30, num.grid.pts.psgz))
for (ff in 1:num.test.models) {
  aero.gz.anom[ff, 1:num.grid.pts.psgz] <- aero.gz.all[ff, 1:num.grid.pts.psgz] - ctl.gz.mn[ptr.CTL[ff], 1:num.grid.pts.psgz]
  aero.gz.anom.1940[ff, 1:num.grid.pts.psgz] <- aero.gz.all.1940[ff, 1:num.grid.pts.psgz] - ctl.gz.mn[ptr.CTL[ff], 1:num.grid.pts.psgz]
  aero.ps.anom[ff, 1:num.grid.pts.psgz] <- 0.01 * (aero.ps.all[ff, 1:num.grid.pts.psgz] - ctl.ps.mn[ptr.CTL[ff], 1:num.grid.pts.psgz])   ###  scale by 0.01 to get units of hPa
  ghg.gz.anom[ff, 1:num.grid.pts.psgz] <- ghg.gz.all[ff, 1:num.grid.pts.psgz] - ctl.gz.mn[ptr.CTL[ff], 1:num.grid.pts.psgz]
  ghg.gz.anom.1940[ff, 1:num.grid.pts.psgz] <- ghg.gz.all.1940[ff, 1:num.grid.pts.psgz] - ctl.gz.mn[ptr.CTL[ff], 1:num.grid.pts.psgz]
  ghg.ps.anom[ff, 1:num.grid.pts.psgz] <- 0.01 * (ghg.ps.all[ff, 1:num.grid.pts.psgz] - ctl.ps.mn[ptr.CTL[ff], 1:num.grid.pts.psgz])   ###  scale by 0.01 to get units of hPa
  nat.gz.anom[ff, 1:num.grid.pts.psgz] <- nat.gz.all[ff, 1:num.grid.pts.psgz] - ctl.gz.mn[ptr.CTL[ff], 1:num.grid.pts.psgz]
  nat.gz.anom.1940[ff, 1:num.grid.pts.psgz] <- nat.gz.all.1940[ff, 1:num.grid.pts.psgz] - ctl.gz.mn[ptr.CTL[ff], 1:num.grid.pts.psgz]
  nat.ps.anom[ff, 1:num.grid.pts.psgz] <- 0.01 * (nat.ps.all[ff, 1:num.grid.pts.psgz] - ctl.ps.mn[ptr.CTL[ff], 1:num.grid.pts.psgz])   ###  scale by 0.01 to get units of hPa
  for (yy in 1:30) { nat.gz.anom.ann[ff, yy, 1:num.grid.pts.psgz] <- nat.gz.all.ann[ff, yy, 1:num.grid.pts.psgz] - ctl.gz.mn[ptr.CTL[ff], 1:num.grid.pts.psgz] }
}

ag.ps.anom<- aero.ps.anom + ghg.ps.anom
ag.gz.anom<- aero.gz.anom + ghg.gz.anom
ag.gz.anom.1940<- aero.gz.anom.1940 + ghg.gz.anom.1940


aero.thick.anom<- aero.gz.anom - 10.0 * aero.ps.anom
ghg.thick.anom<- ghg.gz.anom - 10.0 * ghg.ps.anom
ag.thick.anom<- aero.thick.anom + ghg.thick.anom
nat.thick.anom<- nat.gz.anom - 10.0 * nat.ps.anom


aero.gz.anom.model<- array(0.0, c(num.ctl.models, num.grid.pts.psgz));  aero.ps.anom.model<- aero.gz.anom.model;  ghg.gz.anom.model<- aero.gz.anom.model;  ghg.ps.anom.model<- aero.gz.anom.model
nat.gz.anom.model<- aero.gz.anom.model;  nat.ps.anom.model<- aero.gz.anom.model; aero.gz.anom.model.1940<- aero.gz.anom.model
for (mm in 1:num.ctl.models) {
  ptr.mod<- which(ptr.CTL == mm)
  aero.gz.anom.model[mm, ]<- colMeans(aero.gz.anom[ptr.mod,])
  aero.ps.anom.model[mm, ]<- colMeans(aero.ps.anom[ptr.mod,])
  ghg.gz.anom.model[mm, ]<- colMeans(ghg.gz.anom[ptr.mod,])
  ghg.ps.anom.model[mm, ]<- colMeans(ghg.ps.anom[ptr.mod,])
  nat.gz.anom.model[mm, ]<- colMeans(nat.gz.anom[ptr.mod,])
  nat.ps.anom.model[mm, ]<- colMeans(nat.ps.anom[ptr.mod,])
  aero.gz.anom.model.1940[mm, ]<- colMeans(aero.gz.anom[ptr.mod,]) - colMeans(aero.gz.anom.1940[ptr.mod,])
}

ag.gz.anom.model<- aero.gz.anom.model + ghg.gz.anom.model
ag.ps.anom.model<- aero.ps.anom.model + ghg.ps.anom.model

aero.thick.anom.model<- aero.gz.anom.model - 10.0 * aero.ps.anom.model
ghg.thick.anom.model<- ghg.gz.anom.model - 10.0 * ghg.ps.anom.model
nat.thick.anom.model<- nat.gz.anom.model - 10.0 * nat.ps.anom.model
ag.thick.anom.model<- aero.thick.anom.model + ghg.thick.anom.model



#####  compute the same for ERA5

### sort(sapply(ls(),function(x){object.size(get(x))}))      


ip.fn<- 'D:/***/ps_OctMar_global_grid.csv'
ip.dat<- as.data.frame(fread(ip.fn, header=T, sep=','))
ptr.inc<- which(as.numeric(ip.dat[,2]) > -0.01)
era5.lon<- as.numeric(ip.dat[ptr.inc,1]); era5.lat<- as.numeric(ip.dat[ptr.inc,2]); num.pts.era5<- length(era5.lon); num.yrs.era5<- dim(ip.dat)[2] - 2
era5.ps.anom<- array(0.0, num.pts.era5); era5.ps.anom.9295<- era5.ps.anom; era5.ps.anom.8394<- era5.ps.anom; era5.ps.anom.1322<- era5.ps.anom
for (i in 1:num.pts.era5) {
  temparr<- as.numeric(ip.dat[ptr.inc[i],3:(num.yrs.era5+2)])
  era5.ps.anom[i]<- 0.01 * (mean(temparr[41:60]) - mean(temparr))
  era5.ps.anom.9295[i]<- 0.01 * (mean(temparr[52:55]) - mean(temparr))
  era5.ps.anom.8394[i]<- 0.01 * (mean(temparr[44:55]) - mean(temparr))
  era5.ps.anom.1322[i]<- 0.01 * (mean(temparr[74:83]) - mean(temparr))
}


ip.fn<- 'D:/***/gz500_OctMar_global_grid.csv'
ip.dat<- as.data.frame(fread(ip.fn, header=T, sep=','))
ptr.inc<- which(as.numeric(ip.dat[,2]) > -0.01)
era5.lon<- as.numeric(ip.dat[ptr.inc,1]); era5.lat<- as.numeric(ip.dat[ptr.inc,2]); num.pts.era5<- length(era5.lon); num.yrs.era5<- dim(ip.dat)[2] - 2
era5.gz500.anom<- array(0.0, num.pts.era5); era5.gz500.test<- array(0.0, c(num.pts.era5, num.yrs.era5) ); era5.gz500.anom.9295<- era5.gz500.anom; era5.gz500.anom.4014<- era5.gz500.anom
era5.gz500.anom.8394<- era5.gz500.anom; era5.gz500.anom.1322<- era5.gz500.anom
for (i in 1:num.pts.era5) {
  temparr<- as.numeric(ip.dat[ptr.inc[i],3:(num.yrs.era5+2)])
  era5.gz500.test[i, 1:num.yrs.era5] <- temparr[1:num.yrs.era5]
  era5.gz500.anom[i]<- 0.1 * (mean(temparr[41:60]) - mean(temparr))
  era5.gz500.anom.4014[i]<- 0.1 * (mean(temparr[41:60]) - mean(temparr[1:75]))
  era5.gz500.anom.9295[i]<- 0.1 * (mean(temparr[52:55]) - mean(temparr))
  era5.gz500.anom.8394[i]<- 0.1 * (mean(temparr[44:55]) - mean(temparr))
  era5.gz500.anom.1322[i]<- 0.1 * (mean(temparr[74:83]) - mean(temparr))
}


era5.thick.anom<- era5.gz500.anom - 10.0 * era5.ps.anom
era5.thick.anom.8394<- era5.gz500.anom.8394 - 10.0 * era5.ps.anom.8394
era5.thick.anom.1322<- era5.gz500.anom.1322 - 10.0 * era5.ps.anom.1322


ptr.eu.north.era5<- which(era5.lat >= 60 & era5.lat <= 70 & era5.lon >= -20  & era5.lon <= 10 )
ptr.eu.south.era5<- which(era5.lat >= 35 & era5.lat <= 45 & era5.lon >= -20  & era5.lon <= 10 )
era5.gz500.grad<- array(0.0, num.yrs.era5)
for (y in 1:num.yrs.era5) {
  era5.gz500.grad[y]<- 0.1 * ( mean(era5.gz500.test[ptr.eu.south.era5, y]) - mean(era5.gz500.test[ptr.eu.north.era5, y]) )
}




#######################################################################################################################################

### plot out Figure 3


init.map.plot.2 <- function (fn.xlim = c(-10,30), fn.ylim = c(40,70), fn.xlab = "", fn.ylab = "", fn.cex.axis=1){
  plot(c(0,0), xlim = fn.xlim, ylim = fn.ylim, xaxs = 'i', yaxs = 'i', axes=FALSE,  xlab = "", ylab = "", cex.axis= fn.cex.axis, add = FALSE) }

xlim.usr <- c(-180,180); ylim.usr <- c(0,90)

usr.axis.sz<- 1.75;  usr.leg.size<- 1.37; usr.pt.size<- 0.5
grid.lat.psgz.2<- grid.lat.psgz-0.25;  grid.lat.psgz.3<- grid.lat.psgz+0.2

plotclr <- c('blue4', 'blue1', 'dodgerblue1', 'lightblue3', 'lightblue1', 'green2', 'greenyellow', 'gold1','orange2','red3')
vbreak<- c(-1e8, -30, -24, -18, -12, -6, 0, 6, 12, 18, 1e8)


temparr.20c<- aero.gz.anom 

stderr.aero.gz<- array(0.0, num.grid.pts.psgz);  scale.fac<- 1/sqrt(num.test.models)
for (ii in 1:num.grid.pts.psgz) {  stderr.aero.gz[ii]<- sd(temparr.20c[, ii]) * scale.fac  }

usr.plt.var.mn<- colMeans(temparr.20c)

ptr.signif.aero.gz<- array(0.0, num.grid.pts.psgz)
temparr1<- usr.plt.var.mn + 2*stderr.aero.gz;  px<- which(temparr1 < 0.0); if (length(px) > 0) { ptr.signif.aero.gz[px]<- 1.0 }
temparr2<- usr.plt.var.mn - 2*stderr.aero.gz;  py<- which(temparr2 > 0.0); if (length(py) > 0) { ptr.signif.aero.gz[py]<- 1.0 }



plotclr <- c('blue4', 'blue1', 'dodgerblue1', 'lightblue3', 'lightblue1', 'green2', 'greenyellow', 'gold1','orange2','red3')
vbreak.aero<- c(-1e8, -40, -35, -30, -25, -20, -15, -10, -5, 0, 1e8)

filename<- paste0(opdir.plots, "Figure_3.png")
png(file = filename, width = 2400., height = 1200., res=216, pointsize=10, type = 'cairo')
op<- par( mai=c(0.4, 0.5, 0.2, 0.25), cex=1.0, mfrow=c(1,1))

  init.map.plot.2(fn.xlim = xlim.usr, fn.ylim = ylim.usr, fn.xlab = "", fn.ylab = "", fn.cex.axis=usr.axis.sz)
  colornum <- cut(usr.plt.var.mn, br=vbreak.aero, labels=FALSE); colorcode <- plotclr[colornum]
  points(grid.lon.psgz, grid.lat.psgz, col = colorcode, bg = colorcode, pch=22, cex = usr.pt.size)
  points(grid.lon.psgz, grid.lat.psgz.2, col = colorcode, bg = colorcode, pch=22, cex = usr.pt.size)
  points(grid.lon.psgz, grid.lat.psgz.3, col = colorcode, bg = colorcode, pch=22, cex = usr.pt.size)
  map('world', add=TRUE, lwd=1.05, col='grey5')

  ptr.sig<- which(ptr.signif.aero.gz > 0.5)
  if (length(ptr.sig) > 0) {  usr.pt.size.2<- 0.3*usr.pt.size; points(grid.lon.psgz[ptr.sig], grid.lat.psgz[ptr.sig], col = 'grey18', bg = 'grey18', pch=16, cex = usr.pt.size.2) }

  ax.size = 1.85
  xticks <- seq(-180, 180, by=60);  axis(side=1, at = xticks, labels = xticks, las = 1, tck = -0.015, cex.axis = ax.size, padj=0.05 )
  yticks <- seq(0, 90, by=30);  axis(side=2, at = yticks, labels = yticks, las = 1, tck = -0.015, cex.axis = ax.size, hadj=0.75 )
  
  xx.arr<- unique(grid.lon.psgz);  yy.arr<- unique(grid.lat.psgz); num.usr.rows<- length(xx.arr)
  zz.mat<- matrix(usr.plt.var.mn, nrow=num.usr.rows, byrow=TRUE)
  usrlevs<- c(-40, -35, -30, -25, -20, -15, -10, -5, 0)
  contour(xx.arr, yy.arr, zz.mat, levels = usrlevs, col='white', lty = c(2,2,2,2,2,2,2,2,1), method = "flattest", drawlabels = FALSE, labcex=0.01, add=TRUE)
  text(-26.0, 44.4, '-15', cex=1.75, col='black', pos=4, font=2)
  text(-29.0, 49.44, '-20', cex=1.75, col='black', pos=4, font=2)
  text(-36.0, 54.0, '-25', cex=1.75, col='black', pos=4, font=2)
  text(-41.0, 60.0, '-30', cex=1.75, col='black', pos=4, font=2)

  legend('topleft', legend=leglabs(sprintf('%g',vbreak.aero)), fill=plotclr, cex = usr.leg.size, bty="o")

box(); par(op); dev.off()




##############################



###  end of program
