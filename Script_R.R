library(FactoMineR)
library(factoextra)
library(ggplot2)
library(PerformanceAnalytics)
library(corrplot)
library(lme4)
library(nlme)
library(AER)
library(reshape2)

#### Data manipulation ####
###--kodiak = Geo coordinates
kodiak <- read.table("./data/kodiak.txt")
colnames(kodiak) <- c("Latitude","Longitude")

###--dstns = brabs frequency by size  :
dstns <- read.table("./data/dstns.txt")
colnames(dstns) <- c("Year","Lenght","Females_juvenils","Females_Adult","Males")

####---fleet : data about fisheries :
fleet <- read.table("./data/fleet.txt")
colnames(fleet) <- c("Year","Vessels_register_fishing","Crab_catch","Weight_crab_catch","pot_lifts","Price_crab")
fleet$Year=as.factor(fleet$Year)
fleet$Meankg=fleet$Weight_crab_catch/fleet$Crab_catch
fleet$logprice=log(fleet$Price_crab)

###----fullness : Female frequency and egg laying potential
fullness <- read.table("./data/fullness.txt")
colnames(fullness) <- c("Year","Size","0%_fullness","1-29%_fullnes","30-59%_fullness","60-89%_fullness","90-100%_fullness")

####--survey : The basic survey data
survey <- read.table("./data/survey.txt") 
colnames(survey) = c("annee","district","station","pots","latitude","longitude","pre4","pre3","pre2","pre1","rmale","prmale","juvfemale","adultfemale")
survey$price <- fleet$Price_crab[match(survey$annee,fleet$Year)] 
survey$legal_male <- survey$rmale+survey$prmale 

# New variable : all_male 
survey$all_male <- (survey$pre4+survey$pre3+survey$pre2+survey$pre1+survey$rmale+survey$prmale)
survey$all.female <- survey$juvfemale+survey$adultfemale

#calculation du sex ratio :
survey$sexratio <- survey$all.female/survey$all_male

# 0 --> Nan
survey$sexratio <- replace(survey$sexratio,survey$sexratio==0,NA)
survey$sexratio <- replace(survey$sexratio,survey$sexratio=="NaN",NA)
survey$sexratio <- replace(survey$sexratio,survey$sexratio=="Inf",NA)
survey$district=as.factor(survey$district)
survey$annee=as.factor(survey$annee)

####--eggs: an estimate of the number of eggs per female
eggs <- read.table("./data/eggs.txt",col.names = c("annee","Nboeufs"))
eggs$annee=as.factor(eggs$annee)

###--celsius : ocean temperature at a depth of 100 meters off the Alaskan coast
celsius <- read.table("./data/celsius.txt", col.names = c("annee","mois","temp"))
meanTemp <- aggregate(celsius[,3], list(celsius$annee), mean)
colnames(meanTemp) <- c("annee","temp")
eggs$temp <- meanTemp$temp[match(eggs$annee,meanTemp$annee)]
celsius$annee=as.factor(celsius$annee)
celsius$mois=as.factor(celsius$mois)

###--catch: Commercial catch data for 1960-1982, broken out by district
catch <- read.table("./data/catch.txt",col.names = c("annee","district","Nbcapture","Nbcapturekg"))
catch$annee <- as.factor(catch$annee)
catch$district <- as.factor(catch$district)
catch$Price <- fleet$Price_crab[match(catch$annee,fleet$Year)] 
#ajout du nombre de bateau de pêche commerciale :
catch$vessels <- fleet$Vessels_register_fishing[match(catch$annee,fleet$Year)] 
#ajout du nombre de piège à crabe :
catch$pieges <- fleet$pot_lifts[match(catch$annee,fleet$Year)]

# Price
catch$Pricekg=(catch$Price/0.45)
catch$PoidsCrabs=catch$Nbcapturekg/catch$Nbcapture
catch$Meanprice1crab=catch$PoidsCrabs*catch$Pricekg

####--salinity : ocean salinity at a depth of 100 meters off the Alaskan coast
salinite <- read.table("./data/salinity.txt",col.names = c("annee","mois","salinite"))
meansalinity <- aggregate(salinite[,3], list(salinite$annee), mean)
colnames(meansalinity) <- c("annee","salinity")
eggs$salinity <- meansalinity$salinity[match(eggs$annee,meansalinity$annee)]
salinite$annee=as.factor(salinite$annee)
salinite$mois=as.factor(salinite$mois)

#### descriptive statistics  ####
# Increase in price over time
fleet$Year=as.numeric(fleet$Year)
ggplot(fleet, aes(x = fleet$Year, y= fleet$Price))+
  geom_point(color="darkred", size =3)+
  labs( x = "Annees", y = "Prix en dollard par livre") +
  geom_smooth(model=exp)
fleet$Year=as.factor(fleet$Year)

## Number of captures over time::
ggplot(catch, aes(annee, Nbcapture)) +
  geom_boxplot()+
  labs( x = "Annees", y = "Nombre de crabes capturés", size=5) 
### Number of fishing boats::
ggplot(fleet, aes(x=fleet$Year, y=fleet$Vessels_register_fishing))+
  geom_point()+ 
  labs( x = "Années", y = "Nombre de bateau de pêche enregisté", size=5) 
### Number of traps: :
ggplot(fleet, aes(x=fleet$Year, y=fleet$pot_lifts))+
  geom_point(size=1)+ 
  labs( x = "Années", y = "Nombre de piège à crabe enregisté", size=5) 
### Total weight of catches::
ggplot(fleet, aes(x=fleet$Year, y=fleet$Weight_crab_catch))+
  geom_point(size=1.5)+ 
  labs( x = "Années", y = "Poids total de crabe pêché (en kilogramme)", size=5)
### Average weight of crabs: :
ggplot(catch, aes(annee, PoidsCrabs))+  
  geom_boxplot()+
  labs( x = "Annees", y = "Poids en livre", size=5) 
### Fullness females 100 %
fullness_years <- aggregate(fullness[, 7], list(fullness$Year), sum)   
colnames(fullness_years) <- c("Year",
                              "Female_100%")
ggplot(data=fullness_years,aes(y=fullness_years$`Female_100%`, x=as.factor(fullness_years$Year)))+ 
  geom_point()+
  labs(x = "Années", y = "Nombre de femelle", size=5)

#### Principal Component Analysis (PCA) ####
fleet_stand <- scale(fleet[,2:8])
acp.res2 <- PCA(fleet_stand)
var2 <- get_pca_var(acp.res2)
fviz_eig(acp.res2, addlabels = TRUE, ylim = c(0, 50))
fviz_cos2(acp.res2, choice = "var", axes = 1:2)
acp.res2$eig

fviz_pca_var(acp.res2, col.var = "cos2",axes = 1:2,
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE)  # Évite le chevauchement de texte) 

survey_acp <- cbind.data.frame(survey$latitude,survey$longitude,
                               survey$legal_male,survey$all_male,survey$all.female,survey$sexratio) 
colnames(survey_acp) <- c("latitude","longitude",'legal_mal',"all.male","all.female",
                          "sexratio")
# Standardization of data for PCA
survey_stand <- scale(survey_acp)

acp.res1 <- PCA(survey_stand) 
# missing values in PCA: Replacing NA values with the variable's mean in R. 
var1 <- get_pca_var(acp.res1)
fviz_eig(acp.res1, addlabels = TRUE, ylim = c(0, 50))
fviz_cos2(acp.res1, choice = "var", axes = 1:4)
acp.res1$eig
plot(acp.res1, choix = "var", axes = c(1,2))
corrplot(var1$contrib, is.corr=FALSE) 

# Coloring the representation based on cos²:
fviz_pca_var(acp.res1, col.var = "cos2", axes = c(1,2), 
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE)  # Évite le chevauchement de texte) 

fviz_pca_var(acp.res1, col.var = "cos2", axes = c(1,3), 
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE) 
plot(catch$Meanprice1crab~catch$annee)
hist(catch$Meanprice1crab)

#### Models ####
### Average price of a crab over time
mod_prix1crab1=glm(catch$Meanprice1crab~catch$annee,family = poisson)
plot(mod_prix1crab1)
dispersiontest(mod_prix1crab1)
summary(mod_prix1crab1)
S=summary(mod_prix1crab1)$coefficients

# Simple model:
plot(fleet$Price_crab~fleet$Vessels_register_fishing)
modprix=lm(fleet$Price_crab~fleet$Vessels_register_fishing)
plot(modprix)
# Variance heterogeneity

fleet$logprice=log(fleet$Price_crab)
plot(fleet$logprice~fleet$Vessels_register_fishing+fleet$Meankg)

# Explanatory variables
modprixlog=lm(fleet$logprice~fleet$Vessels_register_fishing+fleet$Crab_catch+fleet$Meankg) 
#### Model used in the report
plot(modprixlog)
modprixlog0=lm(fleet$logprice~1) # Model null
step(modprixlog,scope=~fleet$Vessels_register_fishing+fleet$Crab_catch+fleet$Meankg,direction="both")
summary(modprixlog)
S2=(summary(modprixlog))$coefficients
summary(modprixlog)

 
plot(catch$PoidsCrabs~catch$annee)
modpoids=lm(catch$PoidsCrabs~catch$annee)
plot(modpoids)
summary(modpoids)
S3=summary(modpoids)$coefficients

### Abundance of catchable males
survey$legal_male <- survey$rmale+survey$prmale 
plot(survey$legal_male~survey$annee)

boxplot(survey$all_male~survey$annee)
boxplot(survey$legal_male~survey$annee)

plot(survey$pots~survey$annee)
plot(survey$pots~survey$district)

mod_male=glm(survey$legal_male~survey$annee, family = poisson)
plot(mod_male)

mod_rmale=glm(survey$rmale~survey$annee,family = quasipoisson)
plot(mod_rmale)
summary(mod_rmale)

mod_legalmale=glm(survey$legal_male~survey$annee,family = poisson)
dispersiontest(mod_legalmale)
cor(survey[,4:18])

mod_rmale1=glm(survey$rmale~survey$annee+survey$district,family = quasipoisson)
plot(mod_rmale1)
summary(mod_rmale1)

mod_rmale2=glm(survey$rmale~survey$annee+survey$district+survey$pots,family = quasipoisson)
plot(mod_rmale2)
summary(mod_rmale2)

mod_legalmale0=glm(survey$legal_male~1,family = quasipoisson)
mod_legalmale1=glm(survey$legal_male~survey$annee+survey$district+survey$pots,family = quasipoisson)
survey$annee=relevel(survey$annee, ref = 2) # Changement de reference

# Model used in the report
mod_legalmale=glm(survey$legal_male~survey$annee*survey$district+survey$pots,family = quasipoisson)
plot(mod_legalmale)
summary(mod_legalmale)
anova(mod_legalmale0,mod_legalmale, test = "Chisq")
anova(mod_legalmale1,mod_legalmale, test = "Chisq")
S4=summary(mod_legalmale)$coefficients
pr2=(mod_legalmale$null.deviance-mod_legalmale$deviance)/mod_legalmale$null.deviance
survey$annee=relevel(survey$annee, ref = 2)# Retour ref de base

### Female fertility
meanlegal=aggregate(survey[,16],list(survey$annee),sum)
meansalinity <- aggregate(salinite[,3], list(salinite$annee), mean)
meanfull=aggregate(fullness[,7],list(fullness$Year),sum)
mean29=aggregate(fullness[,4],list(fullness$Year),sum)
survey2=meansalinity
colnames(survey2)=c("Year","salinity")
survey2$temp=meanTemp[1:14,2]
Nbmale=aggregate(survey[,17],list(survey$annee),sum)
Nbfemale=aggregate(survey[,18],list(survey$annee),sum)
survey2=survey2[-1:-3,]
survey2$meansex=Nbfemale[1:11,2]/Nbmale[1:11,2]
survey2$full=meanfull[1:11,2]
survey2$f29=mean29[1:11,2]
survey2$legalmale=meanlegal[1:11,2]

modabond=glm(survey2$legalmale~survey2$full+survey2$f29+survey2$meansex+survey2$temp+survey2$salinity, family=poisson)
dispersiontest(modabond)
modabond=glm(survey2$legalmale~survey2$full+survey2$f29+survey2$meansex+survey2$temp+survey2$salinity, family=quasipoisson)
plot(modabond)
summary(modabond)
S5=summary(modabond)$coefficients
r2=(modabond$null.deviance-modabond$deviance)/modabond$null.deviance

colnames(fullness)=c("Year","Size","0","1-29","30-59","60-89","full")
fullness$Year=as.factor(fullness$Year)
plot(fullness$full~fullness$Year)
mod_full=glm(fullness$full ~fullness$Year, family = poisson)
dispersiontest(mod_full) # SURDISPERSION
mod_full=glm(fullness$full ~fullness$Year, family = quasipoisson)
plot(mod_full)
summary(mod_full)

mod_full2=glm(fullness$full ~fullness$Year+fullness$Size, family = quasipoisson) # Modèle utilisé
plot(mod_full2)
summary(mod_full2)
S6=summary(mod_full2)$coefficients
mod_full60=glm(fullness$`60-89`~fullness$Year+fullness$Size, family = quasipoisson) 
plot(mod_full60)
summary(mod_full60)

meansize=aggregate(fullness[,2],list(fullness$Year),mean)
survey2$size=meansize[1:11,2]

mod_full3=glm(survey2$full~survey2$size+survey2$temp+survey2$salinity,family=quasipoisson)
plot(mod_full3)
summary(mod_full3)
S7=summary(mod_full3)$coefficients

