############################################################################################################



# PREPARATION STEPS



############################################################################################################
# Load the required libraries, and install them if you need to. 
# We will do this with a function that checks if you have a package and installs it if not.
############################################################################################################


install_and_load <- function(package) {
  # Use suppressWarnings to prevent a warning if the package is not installed
  if (!suppressWarnings(require(package, character.only = TRUE))) {
    install.packages(package, dependencies = TRUE)
  }
  library(package, character.only = TRUE)
}


install_and_load("ape")
install_and_load("tidyverse")
install_and_load("ggtree")
install_and_load("caper")
install_and_load("phytools")
install_and_load("nlme")
install_and_load("ggnewscale")
install_and_load("AICcmodavg")


############################################################################################################
# Set the path to your working directory. Mine is shown below, but you will have to change it to yours
############################################################################################################

setwd("~/Dropbox (University of Michigan)/Comparative Workshop Genomics of Migration")

setwd("C:/Users/langebrake/OneDrive - Carl von Ossietzky Universität Oldenburg/Documents/Workshops/MigrationGenomics")

############################################################################################################
# Read in a phylogenetic tree. We are using a tree used by Dufour et al 2024, which can be freely 
# downloaded from here: https://zenodo.org/records/10604034

tree <- read.tree("TreeCharadriiformes_2023_GEBrev_AllSpeciesAdded_consensus.tree")

############################################################################################################



# EXERCISE 1. Plot a tree and rotate its nodes.



############################################################################################################
# The Dufour et al tree has hundreds of species in it. For exercise 1, we are going to select a small number
# of species to visualize
############################################################################################################

tree
# A little bit of info pops up about this tree. 375 tips! That means it shows relatedness among 375 species.
# What species are they?
# This will return a full list of their scientific names:
tree$tip.label


# First, we define the 5 species we will keep for our first plotting exercise.
species_set_1 <- c("Gallinago_delicata", "Gallinago_gallinago", 
                   "Gallinago_paraguaiae", "Gallinago_andina", "Gallinago_nigripennis")

# Now we will define the species we DON'T want to keep, which is the other 370.
# Fortunately we don't need to list them all by hand.

drops1 <- tree$tip.label[!tree$tip.label%in%species_set_1]

# Now we can make a small tree by dropping those 370 extra tips. 
# We'll do that and then make a little plot using base R plotting.

tree1 <- drop.tip(tree, drops1); plot(tree1)

# Base R can make a basic looking tree, but we will use a special package called ggtree to make. 
# Tree plots with the ggplot2 format.
# (If you don't know what that means, don't worry, just follow along).

# To help us follow along when we rotate notes, we can make each species have its own color on the plot.
tip_colors <- c("darkblue","firebrick","darkorchid",
                "seagreen4","darkorange")

labels_colors_df <- data.frame(label = tree1$tip.label, color = tip_colors)



ggtree(tree1) + geom_tiplab(size=6, aes(color=label)) + ggplot2::xlim(0, 10) +
  scale_color_manual(values=setNames(labels_colors_df$color, labels_colors_df$label)) +
  theme(legend.position="none")  
  
# If the species names at the tips are cut off in the plotting window, you can either:
# 1) Make the plot window bigger by hand
# 2) Make the plotting area bigger in the code, like this:
ggtree(tree1) + geom_tiplab(size=6, aes(color=label)) + ggplot2::xlim(0, 20) +
  scale_color_manual(values=setNames(labels_colors_df$color, labels_colors_df$label)) +
  theme(legend.position="none")  
# 3) Make the text smaller, like this:
ggtree(tree1) + geom_tiplab(size=3, aes(color=label)) + ggplot2::xlim(0, 10) +
  scale_color_manual(values=setNames(labels_colors_df$color, labels_colors_df$label)) +
  theme(legend.position="none")  



# ALl of the significant points (tips and nodes) of the tree have numbers, and we will use these numbers 
# to rotate the nodes for this exercise.
# The first set of numbers are the tips. We can see on the plot that there are 5 tips, and we can also
# do this:
Ntip(tree1)

# We can also count the nodes on the tree to see how many there are, or do:
Nnode(tree1)

# The numbers of the nodes are going to be 6-9 (they start counting after the tips).

# We can see how the nodes and tips are numbered like this:
ggtree(tree1) + geom_tiplab(size=6, aes(color=label)) + ggplot2::xlim(0, 10) +
  scale_color_manual(values=setNames(labels_colors_df$color, labels_colors_df$label)) +
  theme(legend.position="none")+ geom_text(aes(label=node), vjust=-1.3, hjust=-0.5)


# A key feature of phylogenetic trees is that their nodes can rotate. The order of the species
# in the vertical orientation doesn't tell you anything about their relationship and their
# evolution. Instead, the branching structure of the tree is what's important.

# The following plots will rotate the tree around its nodes. Notice how the branching structure
# remains the same even while the vertical order of the names on the tree changes.

# Before you run these plots, can you predict how rotating each of the nodes will affect the 
# order of the species names on the right? What will happen when you rotate node 6?
# What about node 9?

ggtree(tree1) %>% rotate(6) + geom_tiplab(size=6, aes(color=label)) + ggplot2::xlim(0, 10) +
  scale_color_manual(values=setNames(labels_colors_df$color, labels_colors_df$label)) +
  theme(legend.position="none")

ggtree(tree1) %>% rotate(7) + geom_tiplab(size=6, aes(color=label)) + ggplot2::xlim(0, 10) +
  scale_color_manual(values=setNames(labels_colors_df$color, labels_colors_df$label)) +
  theme(legend.position="none")  
  
ggtree(tree1) %>% rotate(8) + geom_tiplab(size=6, aes(color=label)) + ggplot2::xlim(0, 10) +
  scale_color_manual(values=setNames(labels_colors_df$color, labels_colors_df$label)) +
  theme(legend.position="none") 

ggtree(tree1) %>% rotate(9) + geom_tiplab(size=6, aes(color=label)) + ggplot2::xlim(0, 10) +
  scale_color_manual(values=setNames(labels_colors_df$color, labels_colors_df$label)) +
  theme(legend.position="none") 



############################################################################################################
# To continue the exercise....
# Using the code above, select a different subset of species to plot.
# Give the species' names different colors.
# You can see ideas for colors here: http://www.stat.columbia.edu/~tzheng/files/Rcolor.pdf 
# Then choose some nodes and try to rotate them.




############################################################################################################


# EXERCISE 2. Ancestral state reconstruction



############################################################################################################

tree <- read.tree("TreeCharadriiformes_2023_GEBrev_AllSpeciesAdded_consensus.tree")

#For ancestral phenotype reconstruction we first need a phenotype!
pheno <- read.csv2("Dufour2019.Phenotypes.csv",header = TRUE)
head(pheno)


#Now we have a huge list of all the birds, and we have to subset our data
#Let`s start with the phenotype data`
PhenoSub <- subset(pheno, Sp.Scien.hbw %in% tree$tip.label)
Ntip(tree)
nrow(PhenoSub)

#hmmm it does not quite match yet! Let us subset the tree according to our dataset

names <- PhenoSub$Sp.Scien.hbw
subtree <- keep.tip(tree, names)
Ntip(subtree)
nrow(PhenoSub)

#yeey! now we can prepare the data for ANCESTRAL STATE RECONSTRUCTION :)
#first, the tree needs to be rooted
is.rooted(subtree)

#Make the tree dichotomous (so each node splits into only two branches)   
tree.dicho <- multi2di(subtree, random=TRUE) 
#Let's plot our new tree
plotTree(tree.dicho,fsize=0.6)

############################################################################################################

#let us bring together the tree and the phenotype data!

##write data to tree, yes we can mae a table out of a tree
t <- as_tibble(tree.dicho)
t


##compile trait phenotype dataset, add any information you want from the data frame
#make sure that you have one column as "label", fitting the tree nomenclature

d <- tibble(label = PhenoSub$Sp.Scien.hbw,
            trait = PhenoSub$strategy_3,
            clade = PhenoSub$Order
)


#combine tree data with phenotype data
td <- full_join(t, d, by = 'label')
head(td)

#now we really start to prepare the reconstruction
#we extract the phenotype info for all tips (row 1-358)
MRP <- td$trait[1:358] #the tree has 358 species and we extract their traits
names(MRP) <- td$label[1:358]

############################################################################################################

#The data is ready! Let`s run the model`
####ace (Maximum Likelihood), model "ER" -> what does it mean?
fit1<-ace(MRP,tree.dicho,model="ER",type="discrete",marginal = TRUE)

##explore and save output, each row shows the likelihood for each phenotype at each node
fit1$lik.anc

TraitRec<-fit1$lik.anc

#isolate the column with the most likely value for each node
node_states = max.col(TraitRec)
node_states<-colnames(TraitRec)[max.col(TraitRec)]

node_states

#insert the node information into our tree-x-data table into the trait column
td[359:715,5] = node_states

############################################################################################################

#plotting!!

(p <- ggtree(tree.dicho, branch.length = "none", size = 1, aes(color=trait)) %<+% td +
    geom_tippoint(aes(colour=trait), cex=2)+
    geom_tiplab(color="black", offset = 0.2, size=1)+
    scale_color_manual(values=c("violet","black","cyan","red"))+
    theme_tree(legend.title = element_text(size = 15),))

#This is a ruff reconstruction of the species ancestral phenotype. If we were very serious, we
#would test different models by specifying different matrices which could explain migration evolution
#instead of the chosen "ER" model and compare the fit. Have a look into the package description
#and find out what other models are out there!

############################################################################################################

#Besides, what do you think about our reconstruction? can we infer the ancestors with high or low 
#certainty? Do you think migration has a high phylogenomic signal?
#actually we do not need to guess, we can just test it!
#we put all the migrants together for simplicity. The tests needs a binary phenotype

d$trait[d$trait == "partial_mig"] = 1
d$trait[d$trait == "strict_mig"] = 1
d$trait[d$trait == "resident"] = 0
d$trait <- as.integer(d$trait)
d<- as.data.frame(d)
Dstat <- phylo.d(d, tree.dicho, label, trait, permut = 10000, rnd.bias=NULL)
Dstat
#there is a plot option!
par(oma=c(3,3,3,3))
plot(Dstat)
#the blue curve represents the model under Brownian process, meaning high phylogenetic signal,
#red comes from random phenotype distribution. Our D is at the black line

#can you find the D statistic? If you consult the help text, what do you conclude for migration 
#evolution in Charadriiformes? What does the plot say?
############################################################################################################




# EXERCISE 3. PGLS



############################################################################################################
# EXERCISE 3.1 -- linear model version
############################################################################################################


# reset the plot parameter from the last exercise
dev.off()

# Here is the dataset on toe dexterity and navigation (FAKE)
trait_data <- read.csv("traits_fake_navigation_with_toes.csv")
# we need to give our dataset row names corresponding to species, so that the tree can match the data
# in our phylogenetic analyses
row.names(trait_data) <- trait_data$Species

# first we will plot the results and see a linear model

ggplot(trait_data, aes(Toe, Nav)) +
  geom_point(size=4) +
  xlab("Toe dexterity") +
  ylab("Navigational precision") + 
  theme_classic(base_size=14) + 
  geom_smooth(method="lm")

linear_model <-  lm(Nav ~ Toe, data = trait_data)
summary(linear_model)

############################################################################################################
# EXERCISE 3.2 -- looking for phylogenetic signal in traits and their residuals, making a PGLS
############################################################################################################



# next we will plot the distribution of the two traits on the tree. This takes place over a series of steps

toedex <- trait_data %>% dplyr::select(Toe)
navprec <- trait_data %>% dplyr::select(Nav)

treeplot <- ggtree(tree)
toeplot <- gheatmap(treeplot, toedex, offset=0, width=0.1, colnames=T,colnames_position = "top",
                    custom_column_labels="", colnames_offset_y = 0.75) +
  scale_fill_viridis_c(option="A", name="toe\ndexterity")
treeplot2 <- toeplot + new_scale_fill()
fullplot <- gheatmap(treeplot2, navprec, offset=6, width=0.1, colnames=T,colnames_position = "top",
                     custom_column_labels="", colnames_offset_y = 0.75) +
  scale_fill_viridis_c(option="D", name="navigational\nprecision")

# show the plot
fullplot 

# Can you see some patterns in the data that look like they correspond to the phylogenetic groups?
# Next we can test this quantitatively
residuals <- resid(linear_model)
phylosig(tree, residuals, method="lambda", test=T)

# What does it mean that the lambda is close to 1? What is the p-value like?


# Next we will do a PGLS to see if the relationship between the traits is driven by the phylogeny

pgls_model <- gls(Nav ~ Toe, correlation = corBrownian(1, phy=tree, form = ~ Species), data = trait_data)
summary(pgls_model)


############################################################################################################
# EXERCISE 3.3 -- looking for phylogenetic signal in traits and their residuals
############################################################################################################



# Use the code below to simulate traits using Brownian Motion along a tree
# Also sample some traits from a normal distribution where each sample is independent
# In each case the two trait sets are simulated independently of each other, but how often do they end up 
# correlated?


trait1.bm <- fastBM(tree)
trait2.bm <- fastBM(tree) 
plot(trait1.bm, trait2.bm, main="Traits simulated with Brownian Motion")

trait1.gauss <- rnorm(375, mean=0, sd=1)
trait2.gauss <- rnorm(375, mean=0, sd=1)
plot(trait1.gauss, trait2.gauss, main="Traits simulated with rnorm")


# run them both a bunch of times and see how they look

# the following code repeats these simulations 100 times each and records whether the traits 
# have a significant p value in a linear model, even though they are not truly related



trait_simulations <- data.frame(bm=rep(NA, 100), gauss=rep(NA, 100))

for(i in 1:100) {
  
  trait1.bm <- fastBM(tree)
  trait2.bm <- fastBM(tree) 
  
  trait_simulations$bm[i] <- summary(lm(trait1.bm~trait2.bm))$coefficients[8]
  
  trait1.gauss <- rnorm(375, mean=0, sd=1)
  trait2.gauss <- rnorm(375, mean=0, sd=1)
  
  trait_simulations$gauss[i] <- summary(lm(trait1.gauss~trait2.gauss))$coefficients[8]
  
}


# take a look at the results. In these tables, the number of cases of "TRUE"
# is the number of times we simulated traits that are correlated just by chance

table(trait_simulations$bm<0.05)
table(trait_simulations$gauss<0.05)

# we can also plot these

hist(trait_simulations$bm, main="Results of 100 simulations\ntraits non-independent (Brownian motion)", xlab="p-values",
     breaks=50)
abline(v=0.05, col="red")
hist(trait_simulations$gauss, main="Results of 100 simulations\ntraits independent ", xlab="p-values",
     breaks=50)
abline(v=0.05, col="red")


# significant p-values occur in the traits that evolve by Brownian Motion tons of times, even though
# we did not make the outcome of one trait dependent on the other!






