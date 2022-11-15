library(ggplot2)

Nosocomial = c(0.378,0.449,0.453,0.424,0.400)
Community = c(0.245,0.305,0.376,0.351,0.320)
Time = c("1997:2000","2001:2004","2005:2008","2009:2012","2013:2016")

df = data.frame(Response = c(Nosocomial,Community),
                Type = c(rep("Nosocomial",5),rep("Community",5)),
                X = rep(Time,2))

ggplot(data = df, aes(x = X, y = Response, color = Type, group = Type)) + 
  geom_line(size = 1.25) + 
  geom_point(size = 4) + 
  xlab("Time period of isolate collection") + 
  ylab ("MRSA Prevalence amongt SA BSI isolates") + 
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=16,face="bold"),
        legend.text = element_text(size = 16),
        legend.title = element_blank())
