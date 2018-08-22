library(oglmx)
library(MASS)
library(plyr)
library(dplyr)
library(ggplot2)
library(ggrepel)
set.seed(153)
data_path = NA #file is available on www.worldvaluessurvey.org/

## validate methods on simulated data
#############
#############

#function to simulate ordinal-scale data 
#assuming two groups and three levels
#given mean and sd of normally distributed latent variables
#assuming reporting cutoffs at 0 and 1
sim_scale_data <- function(n, means, sds, prop_g1 = .5){
    g1 <- rbinom(n,1,prop_g1)
    group <- as.matrix(data.frame("group 1" = g1, "group 2" = 1-g1))
    logsd <- log(sds)
    ystar<- group%*%means + rnorm(n)*exp(group%*%logsd)
    outcome <- rep(0,n)
    outcome[ystar < 0 ] = -1
    outcome[ystar > 1 ] = 1
    simdata <- data.frame(group, outcome)
    simdata
}

# Ordered probit model without accounting for heteroskedasticity
#results.oprob<-oglmx(outcome ~ group, data=simdata, link="probit",
#                     constantMEAN = F, constantSD = F,
#                     delta=0, threshparam = c(0,1))
#summary(results.oprob)

scale_table <- function(data){
    table(data[, c('group.2', 'outcome')])
}

# Ordered probit model accounting for heteroskedasticity
fit_oprobhet <- function(data){
    results.oprobhet<-oglmx(outcome ~ group.2, ~ group.2, data,
                            constantMEAN = T, constantSD = T,
                            threshparam=c(0,1))
    summary(results.oprobhet)
    results.oprobhet
}

#given point estimates for the group means and sds, compute probability 
#that a sample from group 1 is greater than a sample from group 2:
p_greater_fit <- function(modelfit){
    g1_mean = coef(modelfit)[1]
    g2_mean = coef(modelfit)[1] + coef(modelfit)[2]
    g1_sd = exp(coef(modelfit)[3])
    g2_sd = exp(coef(modelfit)[3] + coef(modelfit)[4])
    #p(sample from g1 > sample from g2 | point estimate is correct)
    point_estimate <- pnorm(0, g1_mean - g2_mean, sqrt(g1_sd^2 + g2_sd^2), lower.tail = F)
    
    #using the standard errors for the model fit, you can also compute confidence intervals
    #for the above probability 
    sim = data.frame(mvrnorm(1000, coef(modelfit),  vcov.oglmx(modelfit)))
    sim$g1_mean = sim[, 1]
    sim$g2_mean = sim[, 1] + sim[, 2]
    sim$g1_sd = exp(sim[, 3])
    sim$g2_sd = exp(sim[, 3] + sim[, 4])
    sim <- sim[, c("g1_mean", "g2_mean", "g1_sd", "g2_sd")]
    sampled_probs = pnorm(0, sim$g1_mean - sim$g2_mean, sqrt(sim$g1_sd^2 + sim$g2_sd^2), lower.tail = F)
    conf_bounds <- quantile(sampled_probs, c(.025, .975))
    output = c(conf_bounds[1], point_estimate, conf_bounds[2])
    names(output)[2] = "estimate"
    output
}

p_greater <- function(means, sds){
    pnorm(0, means[1] - means[2], sqrt(sds[1]^2 + sds[2]^2), lower.tail = F)
}


###test calibration of CI
means = c(.4, .44)
sds = c(.7, .9)
n = 5000
in_bounds = 0
num_runs = 100
for (i in 1:num_runs){
    bounds <- p_greater_fit(fit_oprobhet(sim_scale_data(n,means, sds)))[c(1,3)]
    if (i %% 10 == 0){
        print(i)
    } 
    if (p_greater(means, sds) > bounds[1] && p_greater(means, sds) < bounds[2]){
        in_bounds = 1 + in_bounds
    }
}
in_bounds / i

###test fit
means = c(.4, .44)
sds = c(.7, .9)
n = 5000
data <- sim_scale_data(n,means, sds)
scale_table(data)
model = fit_oprobhet(data)
summary(model)
p_greater_fit(model)

## real world data
#############
#############
load(data_path)
wv6_data <- WV6_Data_v_2015_04_18
wv5_data <- WV5_Data_r_v_2015_04_18
x = "8##Albania
12##Algeria
16##American Samoa
20##Andorra
24##Angola
28##Antigua
31##Azerbaijan
32##Argentina
36##Australia
40##Austria
48##Bahrain
50##Bangladesh
51##Armenia
52##Barbados
56##Belgium
60##Bermuda
64##Bhutan
68##Bolivia
70##Bosnia
72##Botswana
76##Brazil
84##Belize
100##Bulgaria
104##Myanmar
108##Burundi
112##Belarus
116##Cambodia
120##Cameroon
124##Canada
144##Sri Lanka
148##Chad
152##Chile
156##China
158##Taiwan
170##Colombia
180##D.R. Congo
184##Cook Islands
188##Costa Rica
191##Croatia
192##Cuba
196##Cyprus (G)
197##Cyprus (T)
203##Czech Rep.
208##Denmark
214##Dominican Rep.
218##Ecuador
222##El Salvador
226##Eq.Guinea
231##Ethiopia
232##Eritrea
233##Estonia
246##Finland
250##France
268##Georgia
270##Gambia
275##Palestine
276##Germany
288##Ghana
292##Gibraltar
300##Greece
320##Guatemala
324##Guinea
328##Guyana
332##Haiti
340##Honduras
344##Hong Kong
348##Hungary
352##Iceland
356##India
360##Indonesia
364##Iran
368##Iraq
372##Ireland
376##Israel
380##Italy
384##Côte d'Ivoire
388##Jamaica
392##Japan
398##Kazakhstan
400##Jordan
404##Kenya
408##North Korea
410##South Korea
414##Kuwait
417##Kyrgyzstan
418##Laos
422##Lebanon
426##Lesotho
428##Latvia
430##Liberia
434##Libya
438##Liechtenstein
440##Lithuania
442##Luxembourg
450##Madagascar
454##Malawi
458##Malaysia
466##Mali
470##Malta
474##Martinique
478##Mauritania
480##Mauritius
484##Mexico
492##Monaco
496##Mongolia
498##Moldova
504##Morocco
508##Mozambique
512##Oman
516##Namibia
524##Nepal
528##Netherlands
554##New Zealand
558##Nicaragua
562##Niger
566##Nigeria
578##Norway
586##Pakistan
591##Panama
598##Papua New Guinea
600##Paraguay
604##Peru
608##Philippines
616##Poland
620##Portugal
624##Guinea-Bissau
626##Timor-Leste
630##Puerto Rico
634##Qatar
642##Romania
643##Russia
646##Rwanda
682##Saudi Arabia
686##Senegal
690##Seychelles
694##Sierra Leone
702##Singapore
703##Slovakia
704##Viet Nam
705##Slovenia
706##Somalia
710##South Africa
716##Zimbabwe
724##Spain
736##Sudan
740##Suriname
752##Sweden
756##Switzerland
760##Syria
762##Tajikistan
764##Thailand
768##Togo
780##Trinidad and Tobago
784##United Arab Emirates
788##Tunisia
792##Turkey
795##Turkmenistan
800##Uganda
804##Ukraine
807##Macedonia
818##Egypt
826##Great Britain
834##Tanzania
840##United States
850##U.S. Virgin Islands
854##Burkina Faso
858##Uruguay
860##Uzbekistan
862##Venezuela
887##Yemen
891##Serbia and Montenegro
894##Zambia
900##West Germany
901##East Germany
902##Tambov
903##Moscow
904##Basque Country
906##Andalusia
907##Galicia
909##North Ireland
910##Valencia
911##Serbia
912##Montenegro
913##Srpska Republic
914##Bosnian Federation
915##Kosovo"
country_codes = sapply(regmatches(x, gregexpr("([[:digit:]]+)##", x))[[1]], function(x){substr(x, 1, nchar(x)-2)}) 
country_names = sapply(regmatches(x, gregexpr("(##[[:alnum:]\\. \\(\\)\\-]+)", x))[[1]], function(x){substr(x, 3, nchar(x))})  #load and parse the country codes
wv5_data$country <- mapvalues(wv5_data$V2, from = country_codes, to = country_names)
wv5_data$happiness <- mapvalues(wv5_data$V10, from = c(1,2,3,4), to = c(4,3,2,1))
happiness_data_wv5 <- wv5_data[wv5_data$happiness %in% c(1,2,3,4), c("country", "happiness")]

results.oprobhet_5<-oglmx(happiness ~ country, ~ country,
                          happiness_data[happiness_data$wave == 5,],
                          constantMEAN = T, constantSD = T,
                          threshparam=c(NA, 1, 2))
summary(results.oprobhet_5)

metrics <- happiness_data[happiness_data$wave == 5,] %>% group_by(country) %>% 
    summarize(
        prop_happy = mean(happiness >= 3),
        prop_very_happy = mean(happiness == 4),
        naive_mean = mean(happiness)
    )

metrics$mean_latent = sapply(metrics$country, FUN = function(x){
    if(x == "Andorra"){
        coef(results.oprobhet_5)[1]
    }
    else{
        coef(results.oprobhet_5)[1] + 
            coef(results.oprobhet_5)[1:(((length(coef(results.oprobhet_5))-3)/2)+1)][[paste("country", x, sep = "")]]
    }
})
metrics$sd_latent = exp(sapply(metrics$country, FUN = function(x){
    if(x == "Andorra"){
        coef(results.oprobhet_5)[1 + ((length(coef(results.oprobhet_5))-3)/2)+1]
    }
    else{
        coef(results.oprobhet_5)[1 + ((length(coef(results.oprobhet_5))-3)/2)+1] + 
            coef(results.oprobhet_5)[(1 + ((length(coef(results.oprobhet_5))-3)/2)+1):(length(coef(results.oprobhet_5))-1)][[paste("country", x, sep = "")]]
    }
}))

#given a fit model, get the "win probability" between two countries
p_happier <- function(country1, country2, modelfit, firstcountry = "Andorra"){
    country1 = paste("country", country1, sep = "")
    country2 = paste("country", country2, sep = "")
    mean_coefs = coef(modelfit)[1:(((length(coef(modelfit))-3)/2)+1)]
    sd_coefs = coef(modelfit)[(1 + ((length(coef(modelfit))-3)/2)+1):(length(coef(modelfit))-1)]
    if (country1 == country2){
        output = c(0.5, 0.5, 0.5)
        names(output) = c("2.5%", "estimate", "97.5%")
        return(output)
    }
    if (country1 == paste("country", firstcountry, sep = "")){
        c1_mean = mean_coefs[1]
        c1_sd = exp(sd_coefs)[1]
    } else{
        c1_mean = mean_coefs[1] + mean_coefs[[country1]]
        c1_sd = exp(sd_coefs[1] + sd_coefs[[country1]])
    }
    if (country2 == paste("country", firstcountry, sep = "")){
        c2_mean = mean_coefs[1]
        c2_sd = exp(sd_coefs)[1]
    } else{
        c2_mean = mean_coefs[1] + mean_coefs[[country2]]
        c2_sd = exp(sd_coefs[1] + sd_coefs[[country2]])
    }
    
    #p(sample from c2 > sample from c1 | point estimate is correct)
    point_estimate <- pnorm(0, c2_mean - c1_mean, sqrt(c1_sd^2 + c2_sd^2), lower.tail = F)
    
    #using the standard errors for the model fit, you can also compute confidence intervals
    #for the above probability 
    sim = data.frame(mvrnorm(5000, coef(modelfit),  vcov.oglmx(modelfit)))
    countries = gsub("[ \\.\\(\\)\\-]", ".", c(country1, country2, firstcountry))
    country1 = countries[1]
    country2 = countries[2]
    firstcountry = countries[3]
    if (country1 == paste("country", firstcountry, sep = "")){
        sim$c1_mean = sim[, 1]
        sim$c1_sd = exp(sim[, 1 + ((length(coef(modelfit))-3)/2)+1 ])
    } else{
        sim$c1_mean = sim[, 1] + sim[, country1]
        sim$c1_sd = exp(sim[, 1 + ((length(coef(modelfit))-3)/2)+1 ] + sim[, paste(country1, "1", sep = '.')])
    }
    if (country2 == paste("country", firstcountry, sep = "")){
        sim$c2_mean = sim[, 1]
        sim$c2_sd = exp(sim[, 1 + ((length(coef(modelfit))-3)/2)+1 ])
    } else{
        sim$c2_mean = sim[, 1] + sim[, country2]
        sim$c2_sd = exp(sim[, 1 + ((length(coef(modelfit))-3)/2)+1 ] + sim[, paste(country2, "1", sep = '.')])
    }
    
    sim <- sim[, c("c1_mean", "c2_mean", "c1_sd", "c2_sd")]
    sampled_probs = pnorm(0, sim$c2_mean - sim$c1_mean, sqrt(sim$c1_sd^2 + sim$c2_sd^2), lower.tail = F)
    conf_bounds <- quantile(sampled_probs, c(.025, .975))
    output = c(conf_bounds[1], point_estimate, conf_bounds[2])
    names(output)[2] = "estimate"
    output
}

#fit the win probabilities for each country pair
for (country1 in metrics$country){
    print(country1)
    for (country2 in metrics$country){
        result = p_happier(country2, country1, results.oprobhet_5)
        metrics[metrics$country == country1, paste(country2, "estimate", sep = "_")] = result[2]
        metrics[metrics$country == country1, paste(country2, "lower", sep = "_")] = result[1]
        metrics[metrics$country == country1, paste(country2, "upper", sep = "_")] = result[3]
    }
}

#plot latent variable distributions by country
#sort by mean
ggplot(metrics, aes(y = reorder(country, mean_latent), xmin = (mean_latent - 1.96*sd_latent), xmax = (mean_latent + 1.96*sd_latent), x = mean_latent )) +
    geom_errorbarh() +
    scale_x_continuous(sec.axis = sec_axis(~ ., name = "Happiness Level",
                                           breaks = c("Not at all                             " = .305,
                                                      "Not very                          " = 1,
                                                      "Rather                                      " = 2, "Very" = 2.395))) +
    labs(x = "Latent Happiness (95% interval)", y = "Country") +
    theme(axis.ticks.x = element_line(colour = c("black", "black", "black", "transparent")))
#sort by lower bound
ggplot(metrics, aes(y = reorder(country, mean_latent - 1.96*sd_latent), xmin = (mean_latent - 1.96*sd_latent), xmax = (mean_latent + 1.96*sd_latent), x = mean_latent )) +
    geom_errorbarh() +
    scale_x_continuous(sec.axis = sec_axis(~ ., name = "Happiness Level",
                                           breaks = c("Not at all                             " = .305,
                                                      "Not very                          " = 1,
                                                      "Rather                                      " = 2, "Very" = 2.395))) +
    labs(x = "Latent Happiness (95% interval)", y = "Country") +
    theme(axis.ticks.x = element_line(colour = c("black", "black", "black", "transparent")))
#sort by upper bound
ggplot(metrics, aes(y = reorder(country, mean_latent + 1.96*sd_latent), xmin = (mean_latent - 1.96*sd_latent), xmax = (mean_latent + 1.96*sd_latent), x = mean_latent )) +
    geom_errorbarh() +
    scale_x_continuous(sec.axis = sec_axis(~ ., name = "Happiness Level",
                                           breaks = c("Not at all                             " = .305,
                                                      "Not very                          " = 1,
                                                      "Rather                                      " = 2, "Very" = 2.395))) +
    labs(x = "Latent Happiness (95% interval)", y = "Country") +
    theme(axis.ticks.x = element_line(colour = c("black", "black", "black", "transparent")))
#sort by variance
ggplot(metrics, aes(y = reorder(country, sd_latent), xmin = (mean_latent - 1.96*sd_latent), xmax = (mean_latent + 1.96*sd_latent), x = mean_latent )) +
    geom_errorbarh() +
    scale_x_continuous(sec.axis = sec_axis(~ ., name = "Happiness Level",
                                           breaks = c("Not at all                             " = .305,
                                                      "Not very                          " = 1,
                                                      "Rather                                      " = 2, "Very" = 2.395))) +
    labs(x = "Latent Happiness (95% interval)", y = "Country") +
    theme(axis.ticks.x = element_line(colour = c("black", "black", "black", "transparent")))

#table of overall metrics
metrics[, 1:6]

#correlation table of overall metrics
cor(metrics[,2:6])

#plot win probability distributions by country for a given reference country
ref_country = "United States"
lower_vname = paste(ref_country, "lower", sep = "_")
upper_vname = paste(ref_country, "upper", sep = "_")
estimate_vname = paste(ref_country, "estimate", sep = "_")
#color-code by happier or unhappier than reference country (at mean and 95% limits of latent distribution)
group = rep("overlapping", dim(metrics)[1])
happier = with(metrics, (get(lower_vname) > .5 & 
                             (mean_latent - sd_latent*1.96) > (mean_latent[country == ref_country] - sd_latent[country == ref_country]*1.96) &
                             (mean_latent + sd_latent*1.96) > (mean_latent[country == ref_country] + sd_latent[country == ref_country]*1.96)))
unhappier = with(metrics, (get(upper_vname) < .5 & 
                               (mean_latent - sd_latent*1.96) < (mean_latent[country == ref_country] - sd_latent[country == ref_country]*1.96) &
                               (mean_latent + sd_latent*1.96) < (mean_latent[country == ref_country] + sd_latent[country == ref_country]*1.96)))
group[happier] = "happier"
group[unhappier] = "unhappier"
group = factor(group, levels =  c("overlapping", "unhappier", "happier") )
if(sum(group=="unhappier") == 0){
    colorvalues = c("#999999", "#66BBFF", "#DD7733")
} else {
    colorvalues=c("#888888", "#DD7733", "#66BBFF")
}
ggplot(metrics, aes(y = reorder(country, get(estimate_vname)), xmin = get(lower_vname), xmax = get(upper_vname), x = get(estimate_vname) )) +
    geom_vline(xintercept = 0.5, color = "DARKGRAY") +
    geom_errorbarh(aes(color = group)) +
    geom_point() +
    scale_x_continuous(limits = c(0,1)) +
    scale_color_manual(values = colorvalues) +
    labs(x = paste("Win Probability over", ref_country, "(95% Confidence Interval)"), y = "Country",
         title = paste("Which Countries Report Being Happier than ", ref_country, "?", sep = ""))

#Distribution plots
m1 = 0
m2 = 10
s1 = 15
s2 = 30
x = seq(-90, 110, .1)
y1 = dnorm(x, m1, s1)
y2 = dnorm(x, m2, s2)
plot(x, y1, type = "l", col = "red", xlab = "x", ylab = "p")
lines(x, y2, col = "blue")

sk = .03
x_right = exp(sk*x)
plot(x_right, y1, type = "l", col = "red", xlab = "x under right-skewed transform", ylab = "p")
lines(x_right, y2, col = "blue")
m1R = exp(sk*m1 + .5*sk*sk*s1*s1)
m2R = exp(sk*m2 + .5*sk*sk*s2*s2)

x_left = -exp(-sk*x)
plot(x_left, y1, type = "l", col = "red", xlab = "x under left-skewed transform", ylab = "p")
lines(x_left, y2, col = "blue")
m1L = -exp(-sk*m1 + .5*sk*sk*s1*s1)
m2L = -exp(-sk*m2 + .5*sk*sk*s2*s2)

#mean x sd plot
ggplot(metrics, aes(x = mean_latent, y = sd_latent, label = country)) +
    geom_point() + 
    geom_text_repel() +
    labs(x = "Mean Latent Happiness", y = "SD Latent Happiness", title = "Two-dimensional model of national happiness")

