#Import flux data
source_1 <- "/Users/alexanderho/Documents/Project_Iso/CSV/Switches_h1_iq.csv"
data_1 = read.csv(source_1, header = FALSE)

source_2 <- "/Users/alexanderho/Documents/Project_Iso/Fly/CSV/Switches_h2_iq.csv"
data_2 = read.csv(source_2, header = FALSE)

#Unpack flux rates
taa.tga.1 <- as.numeric(data_1[2,3])
taa.tga.2 <- as.numeric(data_2[2,3])

tga.taa.1 <- as.numeric(data_1[6,3])
tga.taa.2 <- as.numeric(data_2[6,3])

#Calculate pTGA (see methods for derived formulae)
p.1 <- 1 / (1 + (tga.taa.1 / taa.tga.1))
p.2 <- 1 / (1 + (tga.taa.2 / taa.tga.2))
results <- c(p.1, p.2)

#Add standard deviation bars (as calculated by the python script)
lower <- c(results[1]-0.07639464257336215, results[2]-0.07086761978494859)
upper <- c(results[1]+0.07639464257336215, results[2]+0.07086761978494859)

#Plot
par(mfrow=c(1,1))
barplot(results, xlab = 'Group', ylab='pTGA', ylim=c(0, 1),
        names.arg=c('G+C rich', 'G+C poor'), col=c('gray', 'darkgray'))
arrows(x0=0.7, y0=lower[1], x1=0.7, y1=upper[1], code=3, angle=90, length=0.1)
arrows(x0=1.9, y0=lower[2], x1=1.9, y1=upper[2], code=3, angle=90, length=0.1)
