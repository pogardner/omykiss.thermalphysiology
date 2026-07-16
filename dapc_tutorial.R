## PRACTICE DAPC ###
# tutorial: http://adegenet.r-forge.r-project.org/files/tutorial-dapc.pdf

library(adegenet)
data(dapcIllus)
class(dapcIllus)

names(dapcIllus)

# dapcIllus is a list containing four data sets. We will use a
x <- dapcIllus$a
x

# x is a dataset with 600 individuals simulated under an island model (6 islands) for 30 microsatellite markers
# we use find.clusters to identify clusters, although true clusters are in this case known
pop(x)

# we specify that we want up to k =40 groups (not sure why)
grp <- find.clusters(x, max.n.clust=40)

# the function displays a graph of cumulated variance explained by eigenvalues in teh PCA
# apart from the computational tim, there is no reason for keeping a small number of components; here, we keep all info, so specify to retain
# 200 PCs (there are aroun 110 Pc's)

# next, BIC values are shown for increasing values of k
# this shows that k = 6 is the lowest value, and so we should set k to 6

names(grp)
head(grp$Kstat, 8)
grp$stat
head(grp$grp, 10)
grp$size

# The components are respectively the chosen summary stats (here, BIC) for different values of k
# k, slot Kstat
# selected # of clusters and associated BIC (slot stat)
# group memberships (slot grp)
# group sizes (slot size)

# because we know the actual groups from the dta frame, we can check how well they were retrieved from this procedure
table(pop(x), grp$grp)
table.value(table(pop(x), grp$grp), col.lab=paste("inf", 1:6),
            row.lab=paste("ori", 1:6))

# so this plot shows that the inferred groups agreed with the observed ones! 
# rows correspond to actual groups ("ori") and columns correspond to inferred groups ("inf")
# It's important to note that "k" and the number of clusters are arbitrary. A better question is how many clusters are useful to describe the data?
# BIC plots might not always be clear cut, and that's ok!

