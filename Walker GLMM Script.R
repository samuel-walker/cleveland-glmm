# GLMM for housing vacancy in the Cleveland metropolitan area, 1970-2010
# Samuel Walker (samuel.walker@mail.utoronto.ca)
# Department of Geography
# University of Toronto
# November 2017
# Data from US2010 LTDB (https://s4.ad.brown.edu/projects/diversity/Researcher/Bridging.htm)
# Thanks to Patricio R. Estevez for advice and some R code, see http://tinyurl.com/yd3zsxlk
# Thanks to Ben Bolker for a fantastic FAQ for GLMMs in R: http://bbolker.github.io/mixedmodels-misc/glmmFAQ.html

#########################################
#  Load libraries, load and clean data  #
#########################################

# Set working directory

setwd("C:/Users/Sam/Downloads/R/Output/") 

# Load libraries

library(RCurl)
library(lme4)
library(ggplot2)
library(car)
library(MASS)
library(lmtest)
library(DHARMa)
library(fitdistrplus)
library(Hmisc)
library(texreg)
library(stargazer)
library(htmlTable)
library(reshape2)
library(bbmle)
library(extrafont)

# Load data from Github

mydata = read.csv(text=getURL("https://raw.githubusercontent.com/samuel-walker/cleveland-glmm/master/LTDB_final.csv"), header = TRUE, fileEncoding = "UTF-8")
mydata.full = read.csv("C:/Users/Sam/Downloads/R/LTDB_long.csv", header=TRUE, fileEncoding="UTF-8-BOM")

# Overdispersion parameter estimation function from Ben Bolker: http://bbolker.github.io/mixedmodels-misc/glmmFAQ.html#testing-for-overdispersioncomputing-overdispersion-factor
# Comment additions from http://ase.tufts.edu/gsc/gradresources/guidetomixedmodelsinr/mixed%20model%20guide.html

overdisp_fun <- function(model) {
  ## number of variance parameters in an n-by-n variance-covariance matrix
  vpars <- function(m) {
    nrow(m)*(nrow(m)+1)/2
  }
  # The next two lines calculate the residual degrees of freedom
  model.df <- sum(sapply(VarCorr(model),vpars))+length(fixef(model))
  rdf <- nrow(model.frame(model))-model.df
  # extracts the Pearson residuals
  rp <- residuals(model,type="pearson")
  Pearson.chisq <- sum(rp^2)
  prat <- Pearson.chisq/rdf
  # Generates a p-value. If less than 0.05, the data are overdispersed.
  pval <- pchisq(Pearson.chisq, df=rdf, lower.tail=FALSE)
  c(chisq=Pearson.chisq,ratio=prat,rdf=rdf,p=pval)
}

# Set tract and decade as factors

mydata$TRTID10 <- factor (mydata$TRTID10)
mydata$decade <- factor (mydata$decade)

# Group mean center variables by decade

library(plyr)

head(mydata)

scaled.mydata <- ddply(mydata, c("decade"), transform, P_NONWHT = scale(P_NONWHT), 
                                                       a_hinc = scale(a_hinc))

scaled.mydata <- na.omit(scaled.mydata)

head(scaled.mydata)

# Check missing values

table(is.na(scaled.mydata)) 
sapply(scaled.mydata, function(x) sum(is.na(x)))
table(rowSums(is.na(scaled.mydata)))
table(complete.cases(scaled.mydata))

detach(package:plyr)

# Add ID variable for individual-level random effect later

scaled.mydata$ID <- seq.int(nrow(scaled.mydata))

#########################################
#    Visualizations before modelling    #
#########################################

# Load fonts

font_import(prompt=F)
loadfonts(device = "win") # Note, this can take a long time; remove if you don't care about fonts

# Histogram of dependent

p0 <- ggplot(data=mydata, aes(mydata$R_VAC))+
  geom_histogram(binwidth=250,boundary = 0, closed = "left", fill="gray")+
  labs(x = "Vacant housing units", y = "Count")+
  stat_bin(boundary=0,closed="left",binwidth = 250, geom="text", colour="black", size=2,
           aes(label=..count.., y =..count..*.5)) +
  scale_x_continuous(minor_breaks = NULL, breaks = seq(0,max(mydata$R_VAC), 250))+
  theme_minimal()+
  theme(text=element_text(family="Times New Roman", face="bold", size=12), axis.text.x=element_text(angle = -90, hjust = 0, vjust = 0),
        panel.grid.major.x=element_line())
p0
ggsave("histogram.png")

#position=position_stack(vjust=0.5)

ggplot_build(p0)$data[[1]]$x

#             stat_bin(binwidth=bw, geom="text", colour="white", size=2,
#aes(label=..count.., x=x+bw/2, y=0.5*(..count..))) +

#bw <- diff(range(mydata$R_VAC)) / (2 * IQR(mydata$R_VAC) / length(mydata$R_VAC)^(1/3)) # Freedman-Diaconis rule 
# p0 <- ggplot(data=mydata, aes(mydata$R_VAC))+
#   geom_histogram(binwidth=bw,boundary = 0, closed = "left")+
#   labs(x = "Vacant housing units", y = "Count")+
#   scale_x_continuous(breaks = function(x) unique(floor(seq(0,max(mydata$R_VAC), bw))))+
#   theme_minimal()+
#   theme(text=element_text(family="Times New Roman", face="bold", size=12))
# p0

# Boxplot of vacancy by decade; shows increasing proportion of outliers in later years

p1 <- ggplot(mydata.full, aes(decade, R_VAC/HU, group=decade))+
             geom_point()+
             geom_boxplot(alpha=.5)+
             labs(x = "Decade", y = "Percent vacant housing units")+
             scale_y_continuous(labels = scales::percent)+
             theme_minimal()+
             theme(text=element_text(family="Times New Roman", face="bold", size=12))
p1 
ggsave("boxplot.png")

# Tracing tract vacancy by decade

p2 <- ggplot(mydata.full, aes(decade, R_VAC/HU, group=TRTID10))+
             geom_point(data=mydata.full[mydata.full$D_suburb==1, ], aes(color="b"), size=1)+           
             geom_line(data=mydata.full[mydata.full$D_suburb==1, ], aes(color="b"))+
             geom_point(data=mydata.full[mydata.full$D_suburb==0, ], aes(color="a"), size=1)+
             geom_line(data=mydata.full[mydata.full$D_suburb==0, ], aes(color="a"))+
             labs(x = "Decade", y = "Percent vacant housing units")+
             scale_y_continuous(labels = scales::percent)+
             scale_colour_manual(name = 'Location', labels = c("Cleveland","Other municipalities"), position="left",
                                 values = c("b"="gray","a"="black"))+
             theme_minimal()+
             theme(legend.position=c(.2,.92),legend.title=element_blank(), 
                   text=element_text(family="Times New Roman", face="bold", size=12))
p2
ggsave("tracts_decades.png")

# Percent vacant by decade

library(dplyr)

# Create a group-means data set
gm <- mydata.full %>% 
      group_by_(.dots=c("decade","D_suburb")) %>%
      summarise(
        mean_vac = mean(P_VAC, na.rm=TRUE),
        mean_nonwht = mean(P_NONWHT, na.rm=TRUE),
        mean_ahinc = mean(a_hinc, na.rm=TRUE)
  )
gm

# Plot vacancy
p3 <- ggplot(gm, aes(decade,mean_vac,color=as.logical(D_suburb))) +
             geom_point(data=gm)+
             geom_line(data=gm)+
             labs(x = "Decade", y = "Mean percent vacant housing units by tract")+
             scale_color_manual(labels = c("City of Cleveland", "Other municipalities"), values = c("black","gray"))+
             guides(color=guide_legend(NULL))+
             scale_y_continuous(labels = scales::percent)+
             theme_minimal()+
             theme(legend.position=c(.7,.2),legend.title=element_blank(), 
                   text=element_text(family="Times New Roman", face="bold", size=12))
p3
ggsave("vacancy_by_decade_suburb.png")

# Plot vacancy by decade by place over a multi-page pdf (https://stackoverflow.com/questions/39736655/ggplot2-plots-over-multiple-pages)

pdf("C:/Users/Sam/Downloads/R/Output/facet.pdf", 7, 5)
for (i in seq(1, length(unique(scaled.mydata$place_name)), 6)) {
  print(ggplot(mydata.full[mydata.full$place_name %in% levels(mydata.full$place_name)[i:(i+5)], ], 
               aes(x=decade,y=R_VAC/HU,color=place_name))+
               geom_boxplot(aes(x=decade,y=R_VAC/HU,group=decade))+
               geom_point(aes(x=decade,y=R_VAC/HU))+
               facet_wrap(~place_name)+
               guides(color = "none")+
               labs(x = "Decade", y = "Mean percent vacant housing units by tract")+
               scale_y_continuous(limits=c(0, 0.8), labels = scales::percent)+
               theme_minimal()
               )
}
dev.off()

#########################################
#               Modelling               #
#########################################

# Check conditional mean versus variances as first step in choosing Poisson or negative binomial model

dispersionstats <- scaled.mydata %>%
  group_by(decade) %>%
  summarise(
    means = mean(R_VAC),
    variances = var(R_VAC),
    ratio = variances/means)
dispersionstats

# Variances are much greater than means (see ratio column), indicating that negative binomial is a better choice

# Test Poisson versus negative binomial generalized linear model

modelformula <- formula(R_VAC ~ decade + P_NONWHT * a_hinc + offset(HU_ln))

glm.p <- glm(modelformula, data = scaled.mydata, family = "poisson")   
glm.nb <- glm.nb(modelformula, data = scaled.mydata)

lrtest(glm.p, glm.nb)

# Negative binomial glm is better than Poisson glm

# Test if a mixed model nb (glmm) is more appropriate than a glm nb, just tracts as random

glmmformula <- update(modelformula, . ~ . + (1|TRTID10))

glmm.nb <- glmer.nb(glmmformula, data = scaled.mydata)

anova(glmm.nb, glm.nb) # LRTest not appropriate for GLMM; additionally GLMM must come first in lme4 anova

# Glmm is better

# Test if a nb glmm is better than a Poisson glmm

glmm.p <- glmer(glmmformula, family="poisson", data = scaled.mydata)

anova(glmm.nb, glmm.p)

# Negative binomial glmm is better

# Compare models with new random effects, tracts nested in places

glmmformula <- update(modelformula, . ~ . + (1|place_name/TRTID10))

glmm.nb.place <- glmer.nb(glmmformula, data = scaled.mydata)

anova(glmm.nb, glmm.nb.place)

# Adding place increases log likelihood, but also increases overdispersion.

# Add an individual level random effect to control for overdispersion; makes it a log-normal Poisson

glmmformula <- update(modelformula, . ~ . + (1|ID) + (1|place_name/TRTID10))

glmm.nb.place.id <- glmer.nb(glmmformula, data = scaled.mydata)

anova(glmm.nb.place, glmm.nb.place.id)

# Adding ID reduces AICC by a little bit

# Check overdispersion

overdisp_fun(glmm.nb.place)
overdisp_fun(glmm.nb.place.id)

# Overdispersion not present in model with individual level random effect

# Optimize model fit using allFit, see http://bbolker.github.io/mixedmodels-misc/glmmFAQ.html#troubleshooting

source(system.file("utils", "allFit.R", package="lme4"))
modelfit <- allFit(glmm.nb.place.id)
ss1 <- summary(modelfit)

sink(file="Optimization.txt")
ss1
sink()

sink(file="Optimization.html")
htmlTable(print((modelfit), type="html", useViewer=TRUE))
sink()

summary(modelfit)

# Estimates and log-likelihood all very similar, despite some optimizers throwing warnings. This is ok re: Bolker post! 
# NOTE: need to figure out how to make sure the actual table includes one of the paratmeters that didn't throw a warning. Can't figure that out yet.

# Null nb glmm

glmm.nb.null <- glmer.nb(R_VAC ~ 1 + (1|place_name/TRTID10) + (1|ID) + offset(HU_ln), data = scaled.mydata)

anova(glmm.nb.null, glmm.nb.place.id)

# Present result coefficients as IRRs

fixed <- fixef(glmm.nb.place.id)
confintfixed <- confint(glmm.nb, parm = "beta_", method = "Wald") # Beware: The Wald method is less accurate but much, much faster.
IRR <- exp(cbind(fixed, confintfixed))
IRR

#########################################
#          Model diagnostics            #
#########################################

# Fitted versus residuals

res <- residuals(glmm.nb.place.id)
fit <- fitted(glmm.nb.place.id)

ggplot(data=scaled.mydata,aes(x=fit,y=res))+
       geom_point()+
       geom_hline(yintercept = 0)+
       labs(x="Fitted",y="Residuals")+
       theme_minimal()

# Fitted versus observed

ggplot(data=scaled.mydata,aes(x=fit,y=scaled.mydata$R_VAC))+
  geom_point()+
  geom_smooth(method = "lm", formula = y ~ x)+
  labs(x="Fitted",y="Observed")+
  scale_x_continuous(limits=c(0,1000))+
  scale_y_continuous(limits=c(0,1000))+
  theme_minimal()

# DHARMa

# Generate simulated residuals

simulationOutput <- simulateResiduals(fittedModel = glmm.nb.place.id, n = 1000)

# Plot simulated residuals

png('DHARMa.png', width=1706, height=1278, pointsize = 24)
plotSimulatedResiduals(simulationOutput = simulationOutput)
dev.off()

par(mfrow = c(1,3))
plotResiduals(scaled.mydata$P_NONWHT,  simulationOutput$scaledResiduals)
plotResiduals(scaled.mydata$a_hinc,  simulationOutput$scaledResiduals)
plotResiduals(scaled.mydata$P_NONWHT*scaled.mydata$a_hinc,  simulationOutput$scaledResiduals)
par(mfrow=c(1,1))

# Tests

# K-S test for uniformity of scaled residuals; significant = cannot reject non-uniformity (i.e. evidence of non-uniformity)

sink(file="C:/Users/Sam/Downloads/R/Output/Uniformity.txt")
print(testUniformity(simulationOutput = simulationOutput))
sink()

# Overdispersion cannot do glmer.nb objects, so used Bolker code above instead

# Test for zero inflation (Observed versus expected); significant = cannot reject non-zero inflation (i.e. evidence of zero inflation). 
# Could also be caused by overdipersion, though.

sink(file="C:/Users/Sam/Downloads/R/Output/Zero-inflation.txt")
print(testZeroInflation(simulationOutput))
sink()

# Output simulated residuals to CSV for ArcGIS...

detach(package:dplyr)
library(plyr)

arc <- join(scaled.mydata, simulatedOutput$scaledResiduals, by = "ID", type = "left", match = "all")
head(arc)

arc <- join(arc, simulatedOutput$scaledResidualsNormal, by = "ID", type = "left", match = "all")

detach(package:plyr)

# Pseudo-R2. See http://bbolker.github.io/mixedmodels-misc/glmmFAQ.html. See Xu 2003.

r2.corr.xu <- function(m) {
  1-var(residuals(m))/(var(model.response(model.frame(m))))
}

r2.corr.xu(glmm.p.place.id)

# Different pseudo-R2 formula from Jarret Byrnes on FAQ via http://thread.gmane.org/gmane.comp.lang.r.lme4.devel/684.
# Basically the  correlation between the fitted and observed values

r2.corr.mer <- function(m) {
  lmfit <-  lm(model.response(model.frame(m)) ~ fitted(m))
  summary(lmfit)$r.squared
}

r2.corr.mer(glmm.nb.place.id)

#########################################
#               Outputs                 #
#########################################

# Descriptive statistics

stargazer(scaled.mydata, type = "html", digits=2, 
          out="C:/Users/Sam/Downloads/R/Output/Descriptives.html",
          covariate.labels=c("Vacant housing units","Offset (ln of vacant housing units)","Percent non-white (z-score)","Median household income (z-score)"),
          omit="ID")

stargazer(mydata, type = "html", digits=2, 
          out="C:/Users/Sam/Downloads/R/Output/Descriptives_orig.html",
          covariate.labels=c("Vacant housing units","Offset (ln of vacant housing units)","Percent non-white","Median household income"),
          omit="ID")

# Correlation matricies

sink(file="Correlation_matrix_orig.txt")
rcorr(as.matrix(mydata[,c(4,6,7)]))
sink()

sink(file="Correlation_matrix_scale.txt")
rcorr(as.matrix(scaled.mydata[,c(4,6,7)]))
sink()

# Dispersion statistics as HTML table

sink(file="C:/Users/Sam/Downloads/R/Output/Dispersion_stats.html")
ds <- rapply(object = dispersionstats, f = round, classes = "numeric", how = "replace", digits = 2) 
htmlTable(as.data.frame(ds),header=c("Decade","Means","Variances","Ratio"),rnames=FALSE,align="lccc")
sink()

# AIC results as HTML table (includes AICtable object as text as well for some reason)

AICtable <- AICtab(glm.p, glm.nb, glmm.nb.null, glmm.nb, glmm.p, glmm.nb.place, glmm.nb.place.id,weights=TRUE,delta=TRUE,base=TRUE,
       logLik=TRUE,nobs=TRUE,mnames=list("Poisson generalized linear model","Negative binomial generalized linear model",
       "Negative binomial generalized linear mixed model","Negative binomial generalized linear mixed model with census tract random effects",
       "Poisson generalized linear mixed model","Negative binomial generalized linear model with added census-defined place random effects",
       "Negative binomial generalized linear model with added individual-level random effects"))

sink("AIC.html")
print(htmlTable(print(AICtable),header=c("Log Likelihood","AIC","&#916Log Likelihood",
                "&#916AIC","df","Weight")), type="html",useViewer=TRUE)
sink()

# Overdispersion test results as table

list <- c(glmm.p, glmm.nb.null, glmm.nb, glmm.nb.place, glmm.nb.place.id)
overdisp_rows <- round(mapply(overdisp_fun, list),3)

sink("Overdispersion_table.html")
htmlTable(overdisp_rows, header=c("Model 3","Model 4","Model 5","Model 6","Model 7"),
          rnames=c("Chi-square","Ratio","Residual df","<i>p</i>"))
sink()

# GLMM results as HTML table

htmlreg(list(glmm.nb.null, glmm.nb.place, glmm.nb.place.id), 
        "C:/Users/Sam/Downloads/R/Output/Table.html",
        single.row = TRUE,
        custom.coef.names = c("(Intercept)","1980","1990","2000","2010","Percent non-white population (z-score)",
                              "Median household income (z-score)","% non-white*income"),
        custom.gof.names=c("AIC","BIC","Log Likelihood","Number of observations","Groups: individual-level random effect",
                           "Groups: tracts in census-defined places","Groups: census-defined places",
                           "Variance: individual-level random effect","Variance: tracts in census-defined places","Variance: census-defined places"),
        digits=3,
        caption="")