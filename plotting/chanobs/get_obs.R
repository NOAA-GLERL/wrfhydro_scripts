library(dataRetrieval)
t0 <- '2020-04-01T00:00Z'
tf <- '2020-05-10T00:00Z'

gages <- read.table('gages.txt', colClasses=c('numeric','character'))
OBS <- readNWISuv(siteNumbers=gages[,2], parameterCd='00060', startDate=t0, endDate=tf)
OBS <- OBS[grepl('00:00', OBS$dateTime),] # subset only hourly data
save(OBS, file='OBS.Rdata')

