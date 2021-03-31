#analysis of planktonic foraminiferal d18O to track SO SST front over deglaciation 
#this is an updated version of the method published in Gray, W.R. et al (2020) doi.org/10.1029/2019GL086328 
#it accompanies the manuscript Gray et al (submitted) and has been developed for the Southern Ocean

#mgcv must be installed
require(mgcv)

###
#import compiled d18O data and remove any NAs
dat<-read.csv("SO_d18O_data.csv")
dat<- na.omit(dat)

###
#add unique identifier to each d18O record (i.e. if there are mulitple species analysed onthe same core)
id<-matrix(, nrow = nrow(dat), ncol = 1)

for (i in 1:nrow(dat)){
	id[i]<-paste(dat$core[i],"_",dat$species[i],sep="")
}

dat<-cbind(dat,id)

###
# select the geographic region to over which to perform the analysis
#ALL is ALL, PAC is INDIAN-PACIFIC, ATL is ATLANTIC
#basin<- 'ALL'
basin<- 'PAC' 
#basin<- 'ATL'

#this includes all data south of 60S in the analysis for all basins
dat_pole<- subset(dat, lat_n < -60) #southern cores
dat2<- subset(dat, lat_n > -60) #all data except southern

dat<- if(basin == 'ALL') {dat2} else
      if(basin == 'ATL') {subset(dat2, lon_e > -65 & lon_e < 25)} else
      if(basin == 'PAC') {subset(dat2, lon_e < -65 | lon_e > 25)} else

dat<- rbind(dat,dat_pole)


####
# age settings
age_diff<- 10 #this is the age to difference from, in ka
start_age<- 20 #this is the age to start the analysis, in ka
end_age<- 10 #this is the age to end the analysis, in ka
int<- -0.25 #this is the age interval, in ka
int_age<- seq(start_age,end_age,by=int) #time steps to perform DLat analysis


###
#lat settings
#latitudinal window to minimise differences between time steps (max holocene dd18o/dlat is 44)
lat_s<- -50 #lat south
lat_n<- -40 #lat north
#magnitude and resolution of shift to apply in analysis
delta_lat_min<- -10 #shift min 
delta_lat_max<- 10 #shift max
lat_res <- 0.1 #shift int
kay<- 7 #number of basis functions within lat gam - chosen using gam.check() on LGM and Holocene data


### define shift function, which will be used later to calculate the lat shift between time steps
shift <- function(x, n, invert=FALSE, default=NA){
  stopifnot(length(x)>=n)
  if(n==0){
    return(x)
  }
  n <- ifelse(invert, n*(-1), n)
  if(n<0){
    n <- abs(n)
    forward=FALSE
  }else{
    forward=TRUE
  }
  if(forward){
    return(c(rep(default, n), x[seq_len(length(x)-n)]))
  }
  if(!forward){
    return(c(x[seq_len(length(x)-n)+n], rep(default, n)))
  }
}

###
#import selevel and global temperature data
esl_dat<-read.csv("lambeck2014.csv") #sea level data from lambeck 2014
dt_dat<- read.csv("shakun2012.csv") #global temperature anomaly from shakun 2012
global_mean_sst_DT<- read.csv("pmip_dt.csv")$DT #global LGM-PI SST anomaly from PMIP3/4 models

#calculate total change and interpolate sea level data to time steps
Desl<- max(esl_dat$esl)-min(esl_dat$esl) #total sea level 
slb<-gam(esl~s(age, k=50),data=esl_dat, method='REML') #model sea level with GAM
sl<- predict(slb,newdata=data.frame(age=int_age), se.fit=TRUE) #predict sea level at time steps

#calculate total change and interpolate global temperature data to time steps
Ddt<- max(dt_dat$global)-min(dt_dat$global) #total dt
dtb<-gam(global~s(age, k=50),data=dt_dat, method='REML') #model global temperature data with GAM
dt<- predict(dtb,newdata=data.frame(age=int_age), se.fit=TRUE) #predict temperature data at time steps


###############################
#montecarlo/bootstrap error loop 
mcmc_iterations<- 99 #9999 takes several hours on my machine but is required for stable uncertainites, 99 will give a rough approximation

s1<-matrix(, nrow = length(int_age), ncol = mcmc_iterations) #empty matrix for each iteration of DLat_SST

for (it in 1: mcmc_iterations) {

	#bootstrap
	dat_pole<- subset(dat, lat_n < -60) #take cores south of 60S out of bootstrap
	dat_2<- subset(dat, lat_n > -60)
	dat_2<- dat_2[sample(nrow(dat_2), nrow(dat_2),replace=TRUE), ] #bootstrap sampling of dataset
	dat_i<- rbind(dat_2,dat_pole)
	
	#monte carlo uncertainity for sealevel and global temperature 
	sli<- rnorm(length(sl$fit), sl$fit, sl$se.fit) #uncertainty in the sea level curve
	dti<- rnorm(length(dt$fit), dt$fit, dt$se.fit) #uncertainty in the global temperature curve
	global_ivc_tot_per_mil<- rnorm(1,1,0.05) # whole ocean change in d18Osw from Schrag 2002/Adkins 2001
	global_sst_tot_per_mil<- rnorm(1, mean(global_mean_sst_DT), sd(global_mean_sst_DT))*-0.22 #mean SST change from PMIP3+4/kimoneil1997

	#scale sl/dt over deglaciation
	Dd18Oivc<- global_ivc_tot_per_mil *(sli/Desl) #ice volume contribution at each time step
	Dd18Odt<- global_sst_tot_per_mil *(dti/Ddt) #global mean SST contribution at each time step
	Dd18Oesl<- Dd18Oivc + Dd18Odt #combined contribution at each time step


	#for each core fit d18o with gam as function of age (including random d18o+age error) and predict deglacial d18Oicv-gt at time steps 
	id_sum<-as.vector(unique(dat_i$id))
	r1<-matrix(ncol=1,nrow=length(id_sum)) #empty vector for core lats
	r3<-matrix(ncol=length(int_age),nrow=length(id_sum)) #empty matrix for interpolated d18Oivc-gt data at each timestep

	for (i in 1:length(id_sum)){
		a <- as.character(id_sum[i])
		d<-subset(dat_i, id == a)
		lat<- dat_i$lat[match(a,dat_i$id)]
		r1[i,]<- lat
		#only includes records spanning entire age interval, requiring >10 data points over deglaciation within bootstrap resampling
		if (end_age >= min(d$age, na.rm=TRUE) & start_age <= max(d$age, na.rm=TRUE) & length(d$d18o) >= 11){
			k <- ifelse(length(d$d18o) > 30, 20, 10) #number of basis functions
			#monte carlo error on age and d18o uncertainties
			d$age<- rnorm(length(d$age), d$age, 0.5) #500 yr 1s age uncertainties
			d$d18o<- rnorm(length(d$d18o), d$d18o, 0.04) #0.04 per mil 1s d18o uncertainty 
			#fit gam to record and predict to timestep
			b<-gam(d18o~s(age, k=k),data=d, method='REML')
			int_d18o<-predict(b,newdata=data.frame(age=int_age), se.fit=FALSE)
			r3[i,]<- 	int_d18o + Dd18Oesl } 
			else {r3[i,]<- NA}
		}


	#fit d18o at each time step with gam as function of lat 
	r4<-matrix(ncol=length(int_age),nrow=length(seq(lat_s-delta_lat_max,lat_n-delta_lat_min,by=lat_res))) #empty matrix for lat gams

	for (j in 1:length(int_age)){
		temp<- cbind(r1,r3[,j])
		temp<- na.omit(temp) 
		lat_int <- temp[,1]
		d18o_int <- temp[,2]
		b_int<-gam(d18o_int~s(lat_int, k=kay))
		pdat2<- data.frame(lat_int=seq(lat_s-delta_lat_max,lat_n-delta_lat_min,by=lat_res))
		p_int<-as.numeric(predict(b_int,newdata=pdat2, se.fit=TRUE)$fit)
		r4[,j]<- p_int	
		}


	#delta_lat calculation minimies euclidean distance between gam fit at each timestep and differencing time step
	pdat<- data.frame(lat=seq(lat_s-delta_lat_max,lat_n-delta_lat_min,by=lat_res))
	q_holo<- r4[match(lat_s,pdat$lat):match(lat_n,pdat$lat),match(age_diff,int_age)]
	r5<-matrix(ncol=1,nrow=length(int_age))

	for (j in 1:length(int_age)){
		p_int<- r4[,j]
		delta_lat<- seq(delta_lat_min,delta_lat_max, by=0.1)
		y<-matrix(ncol=1,nrow=length(delta_lat))
		
		for (i in 1:length(delta_lat)){
			nshift<- delta_lat[i]*10
			p_shift <- shift(p_int,nshift)
			q<-p_shift[match(lat_s,pdat$lat):match(lat_n,pdat$lat)]
			y[i,] <- sum((q_holo - q)^2) #calculates euclidean distance at each lat shift
			}
	
		r5[j,]<- -delta_lat[match(min(y),y)] #finds shift with minimum euclidean distance
		}

	s1[,it]<- r5
}

#################################
#bootstrap/montecarlo loop over!



#DLat quantiles
DLat<- apply(s1, 1, quantile, probs=0.5, na.rm=TRUE) #median value
DLat_upr95<- apply(s1, 1, quantile, probs=0.95, na.rm=TRUE) #95% range
DLat_lwr95<- apply(s1, 1, quantile, probs=0.05, na.rm=TRUE) #95% range
DLat_upr68<- apply(s1, 1, quantile, probs=0.68, na.rm=TRUE) #68% range
DLat_lwr68<- apply(s1, 1, quantile, probs=0.32, na.rm=TRUE) #68% range
age<- int_age 

#combine and export results
results<-data.frame(age=age,DLat=DLat,DLat_upr95=DLat_upr95,DLat_lwr95=DLat_lwr95,DLat_upr68=DLat_upr68,DLat_lwr68=DLat_lwr68)
write.table(results,file='DLat_SST_results.csv', sep=',', row.names=FALSE, col.names=TRUE)


#plot results!
dev.new(width=3.5, height=3.5)
par(ps = 10, cex = 1, cex.main = 1); par(mgp=c(1.5,0.5,0)); par(las=1);par(tck=-0.01)
par(mar=c(3.5,3.5,0.5,0.5))

plot(-99, -99, col=adjustcolor("grey99",alpha=0.01),type='o',xlim=c(20,10),ylim=c(6,-1), pch=1, axes=FALSE, xlab='', ylab='')
grid (NULL,NULL, lty = 6, lwd=0.75, col = adjustcolor("cornsilk2", alpha=0.7))

polygon(x=c(age, rev(age)),y=c(DLat_lwr95, rev(DLat_upr95)),col = adjustcolor("grey67",alpha=0.25),border=NA)
polygon(x=c(age, rev(age)),y=c(DLat_lwr68, rev(DLat_upr68)),col = adjustcolor("grey47",alpha=0.25),border=NA)  

points(age, DLat, col=adjustcolor("grey17",alpha=0.9),type='l', pch=1, lwd=1.5)

axis(2, at=seq(6,-1,by=-1), lwd=1)
mtext(expression(paste(Delta*Lat[SST],' (°N)')), side=2, cex.lab=1,line=1.75, las=0, col="black")

axis(1, at=seq(20,10,by=-1), lwd=1)
mtext('Age (ka)', side=1, cex.lab=1,line=1.75, las=0, col="black")

box()

#bye!
