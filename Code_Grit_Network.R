########################################################
# Psychometric Networks of Grit
# Code by: Scarlett Escudero
# Date (version): 2024/November
########################################################

#===============================================================================
# 1. Load dataset and packages
#===============================================================================
library(psych)
library(qgraph)
library(bootnet)
library(mgm)
library(ggplot2)
library(reshape2)
library(responsePatterns)
library(tidyr)
library(smacof)
library(networktools)
library(EGAnet)
library(BDgraph)
library(parallel)
library(dplyr)

grit_data <- read.csv("grit_data.csv")
grit_data <- grit_data[,-1]

#complete_data <- read.csv("data.csv", sep = "")
#complete_data$ID <- seq(1:nrow(complete_data))
#complete_data <- complete_data[, c(3:14,33,35,99)]
#complete_data <- complete_data[complete.cases(complete_data), ]
#grit_data <- complete_data

# Group:
traits <- c("Perseverance","Consistency", "Consistency", "Perseverance","Consistency", "Perseverance", "Consistency", "Consistency","Perseverance","Perseverance","Consistency","Perseverance") # 2 factors according to Duckworth et al. (2007)

#Nodes:
items <- c(
  "Overcome setbacks to conquer an important challenge. ", #GS1
  "New projects sometimes distract me.", #GS2*
  "My interests change from year to year.", #GS3*
  "Setbacks don’t discourage me.", #GS4
  "Obsessed with a project but later lost interest.", #GS5*
  "I am a hard worker.", #GS6
  "Later choose to pursue a different goal.", #GS7*
  "Difficulty maintaining focus.", #GS8*
  "I finish whatever I begin.", #GS9
  "I have achieved a goal that took years of work.", #GS10
  "Interested in new pursuits every few months.", #GS11*
  "I am diligent.") #GS12


set.seed(12345) #Set seed for reproducibility
#===============================================================================
# 2. Descriptives and data cleaning
#===============================================================================

# Aberrant/outlier response patterns
# Mechanistic approach - repetition of response patterns
dataID <- data.frame(ID = c(1:nrow(grit_data)), grit_data)

rppatterns <- rp.patterns(dataID, na.rm=TRUE, id.var = "ID", std.patterns = TRUE, store.data = TRUE)

outliers <- rp.select(rppatterns,percentile=90)@indices[rp.select(rppatterns,percentile=90)@indices$score>0.08,] #cases above percentile 90 AND high scores (>.08)
outliers <- outliers[order(outliers$score),] #428 outliers
as.numeric(rownames(outliers)) #ID that are outliers

for (i in as.numeric(rownames(outliers))) {
  rp.plot(rppatterns,rowname=i,plot=TRUE,
    groups = list(c(2,3,5,7,8,11),c(1,4,6,9,10,12))) #Groups of items 
}
grit_data <- grit_data[!row.names(grit_data) %in% rownames(outliers),] #Remove outliers #dim: 3755 persons

table(grit_data$gender)
mean(grit_data$age, na.rm = TRUE)
sd(grit_data$age, na.rm = TRUE)

grit_data <- grit_data[,-c(13,14,15)]
#write.csv(grit_data,"grit_data.csv")
describe(grit_data)
mardia(grit_data)

# Correlations
datacor <- round(cor_auto(grit_data),2)
melted_datacor <- melt(datacor)
melted_datacor
ggplot(data = melted_datacor, aes(x = Var1, y = Var2, fill = value)) + 
  geom_tile() +
  labs(y = "", x = "")

# Shape of each variable
grit_data_long <- gather(grit_data)
grit_data_long$key <- factor(grit_data_long$key, levels = c("GS1", "GS2", "GS3", "GS4", "GS5", "GS6", "GS7", "GS8", "GS9", "GS10", "GS11", "GS12"))
ggplot(grit_data_long, aes(value)) + 
  geom_bar() + 
  facet_wrap(~key, scales = "free_x") +
  labs(x = "Value", y = "Count") +
  theme(strip.text = element_text(size = 12))

# Search for uninformative nodes
sds <- as.vector(sapply(grit_data, sd, na.rm=T))
psych::describe(sds)
min(sds) #Least informative node
sds < mean(sds)-2.5*sd(sds) #All nodes are informative
sds > mean(sds)+2.5*sd(sds) #All nodes are informative

# Collider bias
results0c <- estimateNetwork(grit_data, default = "cor",
  corMethod = "spearman", 
  threshold = "holm", 
  transform = "rank") # Brute correlatons
results0p <- estimateNetwork(grit_data, default = "pcor",
  corMethod = "spearman",
  threshold = "holm",
  transform = "rank") # Partial correlations

plot0c <- plot(results0c)
plot0p <- plot(results0p)

png("RawCor.png", width = 12, height = 6, res = 600, units = "in")
plot(results0c)
dev.off()
png("PartialCor.png", width = 12, height = 6, res = 600, units = "in")
plot(results0p)
dev.off()

#===============================================================================
# 3. Estimation of the network model
#===============================================================================
# EBIC-glasso
results_EBICglasso <- estimateNetwork(grit_data, 
                                      default = "EBICglasso",
                                      corMethod = "cor_auto",
                                      tuning = 0.5, 
                                      threshold = T)
plot_EBICglasso <- plot(results_EBICglasso)

# ggmModSelect
results_ggmModSelect <- estimateNetwork(grit_data, 
                                        default = "ggmModSelect", 
                                        corMethod = "cor_auto",
                                        tuning = 0.5,
                                        stepwise = TRUE)
plot_ggmModSelect <- plot(results_ggmModSelect)

# huge.npn
results_huge_npn <- estimateNetwork(grit_data, 
                                    default = "huge", 
                                    corMethod = "cor_auto",
                                    tuning = 0.5, 
                                    npn = TRUE, 
                                    criterion = "ebic")
plot_huge_npn <- plot(results_huge_npn)


# TMFG
results_TMFG <- estimateNetwork(grit_data, 
                                default = "TMFG",
                                corMethod = "cor_auto",
                                graphType = "pcor", 
                                tuning = 0.5)
plot_TMFG <- plot(results_TMFG)

# mgm.EBIC
results_mgmEBIC <- estimateNetwork(grit_data, 
                                   default = "mgm",
                                   tuning = 0.25,
                                   corMethod = "cor_auto",
                                   criterio = "EBIC",
                                   rule = "OR")
plot_mgmEBIC <- plot(results_mgmEBIC)

# mgm.CV
results_mgmCV <- estimateNetwork(grit_data, 
                                 default = "mgm",
                                 tuning = 0.25,
                                 corMethod = "cor_auto",
                                 criterio = "CV",
                                 nfolds = 10,
                                 rule = "OR")
plot_mgmCV <- plot(results_mgmCV)

# Sparsity indices
weights_EBICglasso <- results_EBICglasso$graph[lower.tri(results_EBICglasso$graph)]
weights_ggmModSelect <- results_ggmModSelect$graph[lower.tri(results_ggmModSelect$graph)]
weights_huge_npn <- results_huge_npn$graph[lower.tri(results_huge_npn$graph)]
weights_TMFG <- results_TMFG$graph[lower.tri(results_TMFG$graph)]
weights_mgmEBIC <- results_mgmEBIC$graph[lower.tri(results_mgmEBIC$graph)]
weights_mgmCV <- results_mgmCV$graph[lower.tri(results_mgmCV$graph)]

sparsities <- cbind(c("EBICglasso","ggmModSelect","huge.npn","TMFG","mgm.EBIC","mgm.CV"),
                    c(length(weights_EBICglasso[weights_EBICglasso==0])/length(weights_EBICglasso),
                      length(weights_ggmModSelect[weights_ggmModSelect==0])/length(weights_ggmModSelect),
                      length(weights_huge_npn[weights_huge_npn==0])/length(weights_huge_npn),
                      length(weights_TMFG[weights_TMFG==0])/length(weights_TMFG),
                      length(weights_mgmEBIC[weights_mgmEBIC==0])/length(weights_mgmEBIC),
                      length(weights_mgmCV[weights_mgmCV==0])/length(weights_mgmCV)))
tibble(sparsities)

averageLayout <- averageLayout(plot_EBICglasso, plot_ggmModSelect, plot_huge_npn, plot_TMFG, plot_mgmEBIC, plot_mgmCV)

png("SixAlgorithms.png", width = 15, height = 8, res = 600, units = "in")
par(mfrow=c(2,3), cex=1.3)
plot(results_EBICglasso,layout = averageLayout, title = "EBICglasso")
plot(results_ggmModSelect,layout = averageLayout, title = "ggmModSelect")
plot(results_huge_npn,layout = averageLayout, title = "HUGE-NPN")
plot(results_TMFG,layout = averageLayout, title = "TMFG")
plot(results_mgmEBIC,layout = averageLayout, title = "mgm-EBIC")
plot(results_mgmCV,layout = averageLayout, title = "mgm-CV")
par(mfrow=c(1,1))
dev.off()

summary(results_TMFG)
#===============================================================================
# 4. Description of the network model
#===============================================================================
results_TMFG$graph #weights

# Fruchterman-Reingold plot
png("Network_Fruchterman-ReingoldPlot.png", width = 12, height = 8, res = 600, units = "in")
plot(results_TMFG, color = c("#fee090","#e0f3f8"))
dev.off()

# MDS plot
dissimilarity_net1 <- sim2diss(cor_auto(grit_data))
net1_MDS <- mds(dissimilarity_net1)
head(round(net1_MDS$conf, 2))
# Check which method of MDS is best for estimating the plot
net1_MDS_ordinal <- mds(dissimilarity_net1, type = "ordinal")
net1_MDS_ratio <- mds(dissimilarity_net1, type = "ratio")
net1_MDS_interval <- mds(dissimilarity_net1, type = "interval")
net1_MDS_mspline <- mds(dissimilarity_net1, type = "mspline")

par(mfrow=c(2,2))
plot(net1_MDS_ordinal, plot.type = "Shepard", main="Ordinal", pch = 19, cex = 1); text(0.7,1, col="red", paste("Stress =",round(net1_MDS_ordinal$stress,2)))
plot(net1_MDS_ratio, plot.type = "Shepard", main="Ratio", pch = 19, cex = 1); text(0.7,1, col="red",paste("Stress =",round(net1_MDS_ratio$stress,2)))
plot(net1_MDS_interval, plot.type = "Shepard", main="Interval", pch = 19, cex = 1); text(0.7,1,col="red",paste("Stress =",round(net1_MDS_interval$stress,2)))
plot(net1_MDS_mspline, plot.type = "Shepard", main="Spline", pch = 19, cex = 1); text(0.7,1, col="red", paste("Stress =",round(net1_MDS_mspline$stress,2)))
par(mfrow=c(1,1)) #Best is spline


png("Network_MDSPlot.png", width = 15, height = 8, res = 600, units = "in")
MDSnet <- MDSnet(results_TMFG$graph,
                 type = c("mspline"),
                 MDSadj = cor_auto(grit_data),
                 #groups=gr2,
                 color=c("#fee090","#e0f3f8"),
                 edge.color = c("red", "blue"),
                 details = FALSE,
                 legend = FALSE,
                 nodeNames = colnames(grit_data),
                 vsize = 5)
dev.off()


# Detection of redundant (very collinear or locally dependent) nodes using UVA
UVA(grit_data,
    network = results_TMFG$graph,
    reduce = TRUE,
    uva.method = "MBR",
    reduce.method = "mean") #usar "latent" si >3

# Centrality
centralityTable(results_TMFG, standardized = TRUE)
centralityPlot(results_TMFG, include = "all", theme_bw = FALSE)


# Bridge centrality
bridgecentrality <- bridge(results_TMFG$graph, traits, useCommunities="all")
plot(bridgecentrality, zscore = TRUE, plotNA = TRUE) 

# R2 and predictive accuracy
library(mgm)
type <- rep('g', results_TMFG$nNode)
level <- rep(1,results_TMFG$nNode)
datamgm <- as.matrix(grit_data[complete.cases(grit_data),])
fit1 <- mgm(datamgm,
          type = type, 
          level = level)

pred1 <- predict(fit1, datamgm)
pred1$error$R2
mean(pred1$error$R2)
sd(pred1$error$R2)
graph1 <- plot(results_TMFG, pie = pred1$error$R2)

# EGA
datacor <- cor_auto(grit_data)
set.seed(42)
ega.walktrap <- EGA(grit_data, model = "glasso", plot.EGA = TRUE); ega.walktrap
ega.walktrapTMFG <- EGA(grit_data, model = "TMFG", plot.EGA = TRUE); ega.walktrapTMFG

for(i in 3:20) {
  sgc <- EGA(datacor, model = "glasso",n = nrow(grit_data), steps = i, plot.EGA = FALSE)
  print(sgc$dim.variables$dimension)
  print(table(sgc$dim.variables$dimension))
}

vn.entropy(data = datacor, structure = ega.walktrap$wc) # Entropy index

png("EGA_plot.png", width = 8, height = 6, res = 600, units = "in")
EGA(grit_data, model = "glasso", plot.EGA = TRUE); ega.walktrap
dev.off()

plot(results_TMFG,
     groups = traits,
     layout = "spring", 
     legend.cex = 0.45,
     vsize = 4,
     border.width = 2, 
     border.color = '#555555',
     minimum = 0.05,
     color = c("#fee090", "#e0f3f8"),
     labels = colnames(data),
     nodeNames = items,
     legend = TRUE)

#===============================================================================
# 4. Replicability simulation study
#===============================================================================

### Pre hoc power analysis ###
# No clusters in the network
Sim <- netSimulator(
  input = function()bootnet::genGGM(12, p = 0.5, graph = "random"),  
  dataGenerator = ggmGenerator(),
  nCases = c(500,2500,3000,5000,10000),
  nCores = detectCores()-2,
  nReps = 500,
  default = "TMFG")

data_sim <- as.data.frame(Sim)
data_long <- data_sim %>%
  pivot_longer(
    cols = c(sensitivity, specificity, correlation),
    names_to = "Metric",
    values_to = "Value")
Sim_plot <- ggplot(data_long, aes(x = factor(nCases), y = Value)) +
  geom_boxplot(outlier.size = 0.5) +
  facet_wrap(~ Metric, scales = "free_y", 
             labeller = labeller(Metric = c(
               correlation = "Correlation",
               sensitivity = "Sensitivity",
               specificity = "Specificity"))) +
  ylim(0,1) +
  labs(x = "Sample size",
    y = "") +
  theme(strip.text = element_text(size = 14),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"); Sim_plot

png("Sim_Plot.png", width = 10, height = 6, res = 600, units = "in")
Sim_plot
dev.off()

data_long2 <- data_sim %>%
  pivot_longer(
    cols = c(strength, closeness, betweenness, ExpectedInfluence),
    names_to = "Metric",
    values_to = "Value")
Sim_plot2 <- ggplot(data_long2, aes(x = factor(nCases), y = Value)) +
  geom_boxplot(outlier.size = 0.5) +
  facet_wrap(~ Metric, scales = "free_y", ncol=4,
             labeller = labeller(Metric = c(
               strength = "Strength",
               closeness = "Closeness",
               betweenness = "Betweenness",
               ExpectedInfluence = "Expected Influence"))) +
  ylim(0,1) +
  labs(x = "Sample size",
       y = "") +
  theme(strip.text = element_text(size = 14),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none"); Sim_plot2
png("Sim_Plot2.png", width = 12, height = 6, res = 600, units = "in")
Sim_plot2
dev.off()

### Pre hoc replicability analysis ###
SimRep1 <- replicationSimulator(input = function()bootnet::genGGM(Nvar = 12, p = 0.5, 
                                                                  graph = "random"),  
                                dataGenerator = ggmGenerator(),
                                nCases = c(500,2500,3000,5000,10000),
                                nCores = detectCores() - 2,
                                nReps = 500,
                                default = "TMFG")

data_SimRep1 <- as.data.frame(SimRep1)
data_long_SimRep1 <- data_SimRep1 %>%
  pivot_longer(
    cols = c(correlation, jaccard, replicatedEdges, replicatedZeroes),
    names_to = "Metric",
    values_to = "Value")

SimRep1_plot <- ggplot(data_long_SimRep1, aes(x = factor(nCases), y = Value)) +
  geom_boxplot(outlier.size = 0.5) +
  facet_wrap(~ Metric, scales = "free_y", ncol = 4, 
             labeller = labeller(Metric = c(
               correlation = "Correlation",
               jaccard = "Jaccard Index",
               replicatedEdges = "Replicated Edges",
               replicatedZeroes = "Replicated Zeros"))) +
  ylim(0,1) +
  labs(x = "Sample size",
       y = "") +
  theme(strip.text = element_text(size = 14),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none"); SimRep1_plot

png("SimRep1_Plot.png", width = 10, height = 6, res = 600, units = "in")
SimRep1_plot
dev.off()

data_long_SimRep2 <- data_SimRep1 %>%
  pivot_longer(
    cols = c(strength, closeness, betweenness),
    names_to = "Metric",
    values_to = "Value")
SimRep1_plot2 <- ggplot(data_long_SimRep2, aes(x = factor(nCases), y = Value)) +
  geom_boxplot(outlier.size = 0.5) +
  facet_wrap(~ Metric, scales = "free_y", ncol = 4,
             labeller = labeller(Metric = c(
               strength = "Strength",
               closeness = "Closeness",
               betweenness = "Betweenness"))) +
  ylim(0,1) +
  labs(x = "Sample size",
       y = "") +
  theme(strip.text = element_text(size = 14),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none"); SimRep1_plot2
png("SimRep1_Plot2.png", width = 12, height = 6, res = 600, units = "in")
SimRep1_plot2
dev.off()

### Case-drop Bootstrap Correlation Stability ###
boot2 <- bootnet(results_TMFG, 
                 nCores = detectCores() - 2,
                 nBoots = 1000, 
                 type ='case',
                 statistics = c("edge","strength","expectedInfluence","Betweenness","Closeness"))

png("Bootstrap_Edge.png", width = 6, height = 4, res = 600, units = "in")
plot(boot2, statistics = "edge")
dev.off()

png("Bootstrap2.png", width = 6, height = 4, res = 600, units = "in")
p2 <- plot(boot2, statistics = c("strength","expectedInfluence","betweenness","closeness"));p2
dev.off()

# Compute CS-coefficients: above .50 ideal
corStability(boot2, cor = 0.7, statistics = "edge", verbose = TRUE) 
corStability(boot2, cor = 0.7, statistics = c("strength","expectedInfluence","betweenness","closeness"), verbose = TRUE) 

