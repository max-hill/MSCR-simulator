#!/usr/bin/env Rscript
##
## make-plots.R --- use the R package ggplot2 to create plots from simulated data
##
## Author: Max Hill
## (Last updated 2021-04-30)





##______________________________________________________________________________
##
## Setup -- Install and load necessary R packages.
##______________________________________________________________________________


packages <- c("ggplot2")

installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages]) # Install packages
}

invisible(lapply(packages, library, character.only = TRUE)) # Packages loading


setwd('../data/') # Set working directory to 'data/' 





##______________________________________________________________________________
##
## Part 1 -- Effect of recombination when all other variables are fixed and
##           internal branch length is small.
##______________________________________________________________________________

## Simulation parameters:
## (defparameter *τ_ab-values* '(1))
## (defparameter *f-values* '(.01))
## (defparameter *ρ_a-values* '(0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20))
## (defparameter *ρ_b-values* '(0))
## (defparameter *ρ_c-values* '(0))
## (defparameter *ρ_ab-values* '(0))
## (defparameter *ρ_abc-values* '(0))
## (defparameter *θ-values* '(.1)) 
## (defparameter *N* 15000)
## (defparameter *L* 500)
## (defparameter *τ_max* 999999)

## jc-expected
df = read.csv("jc-expected-N15000-L500.csv")
plot1 = ggplot(df, aes(x=ρ_a, y= P.bc - P.ab)) +
    geom_point() + 
    labs(title="Recombination effect on most likely inferred gene topology under JC-expected", x="recombination rate in population A",y="Pr[A(BC)]-Pr[(AB)C]")
ggsave("plot1-jce.jpeg",path="../analysis/")

## jc-sequence
df = read.csv("jc-sequence-N15000-L500.csv")
plot1jce = ggplot(df, aes(x=ρ_a, y= P.bc - P.ab)) +
    geom_point() + 
    labs(title="Recombination effect on most likely inferred gene topology under JC-sequence",x="recombination rate in population A",y="Pr[A(BC)]-Pr[(AB)C]")
ggsave("plot1-jcs.jpeg",path="../analysis/")

## ml-sequence
df = read.csv("ml-sequence-N15000-L500.csv")
plot1jce = ggplot(df, aes(x=ρ_a, y= P.bc - P.ab)) +
    geom_point() + 
    labs(title="Recombination effect on most likely inferred gene topology under ML-sequence",x="recombination rate in population A",y="Pr[A(BC)]-Pr[(AB)C]")
ggsave("plot1-mls.jpeg",path="../analysis/")


##______________________________________________________________________________
##
## Part 2 -- Relative effects of internal branch length and recombination rate
##           on correct topological inference for gene triplets.
##______________________________________________________________________________

## Simulation parameters:
## (defparameter *τ_ab-values* '(1))
## (defparameter *f-values* '(.001 .01 .1))
## (defparameter *ρ_a-values* '(0 5 10))
## (defparameter *ρ_b-values* '(0))
## (defparameter *ρ_c-values* '(0))
## (defparameter *ρ_ab-values* '(0 5 10))
## (defparameter *ρ_abc-values* '(0)) 
## (defparameter *θ-values* '(.001 .01 .05 .1 .2))
## (defparameter *N* 10000)
## (defparameter *L* 500)
## (defparameter *τ_max* 999999)

## jc-expected
df = read.csv("jc-expected-N10000-L500.csv")
plot2 = ggplot(df, aes(x=ρ_a + (τ_abc-τ_ab)*ρ_ab/2 , y= P.bc - P.ab , color=τ_abc-τ_ab)) +
    geom_point() + 
    labs(title="Recomb., internal branch, and gene topology probabilities under JC-expected", x="ρ_a + (τ_abc-τ_ab)*ρ_ab/2",y="Pr[A(BC)]-Pr[(AB)C]")
ggsave("plot2-jce.jpeg",path="../analysis/")

## jc-sequence
df = read.csv("jc-sequence-N10000-L500.csv")
plot2 = ggplot(df, aes(x=ρ_a + (τ_abc-τ_ab)*ρ_ab/2 , y= P.bc - P.ab , color=τ_abc-τ_ab)) +
    geom_point() + 
    labs(title="Recomb., internal branch, and gene topology probabilities under JC-sequence", x="ρ_a + (τ_abc-τ_ab)*ρ_ab/2",y="Pr[A(BC)]-Pr[(AB)C]")
ggsave("plot2-jcs.jpeg",path="../analysis/")

## ml-sequence
df = read.csv("ml-sequence-N10000-L500.csv")
plot2 = ggplot(df, aes(x=ρ_a + (τ_abc-τ_ab)*ρ_ab/2 , y= P.bc - P.ab , color=τ_abc-τ_ab)) +
    geom_point() + 
    labs(title="Recomb., internal branch, and gene topology probabilities under ML-sequence", x="ρ_a + (τ_abc-τ_ab)*ρ_ab/2",y="Pr[A(BC)]-Pr[(AB)C]")
ggsave("plot2-mls.jpeg",path="../analysis/")

## Note for future: redo this simulation with more values for f, more values for recombination rate in A, no recombination in AB, and only one mutation rate


##______________________________________________________________________________
##
## Part 3 -- Identifying the inconsistency zones with gradient illustrations.
##______________________________________________________________________________

## Simulation parameters:
## (defparameter *τ_ab-values* '(1))
## (defparameter *f-values* '(.01 .02 .03 .04 .05 .06 .07 .08 .09 .1 .11 .12 .13 .14 .15 .16 .17 .18 .19 .2))
## (defparameter *ρ_a-values* '(0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20))
## (defparameter *ρ_b-values* '(0))
## (defparameter *ρ_c-values* '(0))
## (defparameter *ρ_ab-values* '(0))
## (defparameter *ρ_abc-values* '(0)) 
## (defparameter *θ-values* '(.1 .01)) 
## (defparameter *N* 10000) 
## (defparameter *L* 50) 
## (defparameter *τ_max* 999999)

## ML-sequence
df = subset(read.csv("ml-sequence-N10000-L50.csv"), θ==0.01)
plot3 = ggplot(df, aes(x = ρ_a, y = τ_abc-τ_ab, color = P.bc - P.ab)) +
    geom_point() + 
    labs(title="Gradient illustration of inconsistency zone under ML-sequence with θ=0.01", x="ρ_a",y="f=τ_abc-τ_ab")
ggsave("plot3-mls-th0.01.jpeg",path="../analysis/")

df = subset(read.csv("ml-sequence-N10000-L50.csv"), θ==0.1)
plot3 = ggplot(df, aes(x = ρ_a, y = τ_abc-τ_ab, color = P.bc - P.ab)) +
    geom_point() + 
    labs(title="Gradient illustration of inconsistency zone under ML-sequence with θ=0.1", x="ρ_a",y="f=τ_abc-τ_ab")
ggsave("plot3-mls-th0.1.jpeg",path="../analysis/")

## JC-expected
df = subset(read.csv("jc-expected-N10000-L50.csv"), θ==0.01)
plot3 = ggplot(df, aes(x = ρ_a, y = τ_abc-τ_ab, color = P.bc - P.ab)) +
    geom_point() + 
    labs(title="Gradient illustration of inconsistency zone under JC-expected with θ=0.01", x="ρ_a",y="f=τ_abc-τ_ab")
ggsave("plot3-jce-th0.01.jpeg",path="../analysis/")

df = subset(read.csv("jc-expected-N10000-L50.csv"), θ==0.1)
plot3 = ggplot(df, aes(x = ρ_a, y = τ_abc-τ_ab, color = P.bc - P.ab)) +
    geom_point() + 
    labs(title="Gradient illustration of inconsistency zone under JC-expected with θ=0.1",x="ρ_a",y="f=τ_abc-τ_ab")
ggsave("plot3-jce-th0.1.jpeg",path="../analysis/")

## JC-sequence
df = subset(read.csv("jc-sequence-N10000-L50.csv"), θ==0.01)
plot3 = ggplot(df, aes(x = ρ_a, y = τ_abc-τ_ab, color = P.bc - P.ab)) +
    geom_point() + 
    labs(title="Gradient illustration of inconsistency zone under JC-sequence with θ=0.01", x="ρ_a",y="f=τ_abc-τ_ab")
ggsave("plot3-jcs-th0.01.jpeg",path="../analysis/")

df = subset(read.csv("jc-sequence-N10000-L50.csv"), θ==0.1)
plot3 = ggplot(df, aes(x = ρ_a, y = τ_abc-τ_ab, color = P.bc - P.ab)) +
    geom_point() + 
    labs(title="Gradient illustration of inconsistency zone under JC-sequence with θ=0.1", x="ρ_a",y="f=τ_abc-τ_ab")
ggsave("plot3-jcs-th0.1.jpeg",path="../analysis/")


##______________________________________________________________________________
##
## Part 4 -- Discrete illustrationss of inconsistency zone for inference of rooted
##           triplet topology.
##______________________________________________________________________________

## parameters same as for plot 3

# JC-expected
df = subset(read.csv("jc-expected-N10000-L50.csv"), θ==0.1)
p <- df$P.bc>df$P.ab
q <- df$P.ac>df$P.ab
inconsistency_zone <- p | q # the vertical bar | is a vectorized logical 'or'
plot3 = ggplot(df, aes(x = ρ_a, y = τ_abc-τ_ab, color = inconsistency_zone)) +
    geom_point() + 
    labs(title="Approximate inconsistency zone under JC-expected for θ=0.1",x="ρ_a",y="f=τ_abc-τ_ab")
ggsave("plot4-jce-th0.1.jpeg",path="../analysis/")

df = subset(read.csv("jc-expected-N10000-L50.csv"), θ==0.01)
p <- df$P.bc>df$P.ab
q <- df$P.ac>df$P.ab
inconsistency_zone <- p | q
plot3 = ggplot(df, aes(x = ρ_a, y = τ_abc-τ_ab, color = inconsistency_zone)) +
    geom_point() + 
    labs(title="Approximate inconsistency zone under JC-expected for θ=0.01",x="ρ_a",y="f=τ_abc-τ_ab")
ggsave("plot4-jce-th0.01.jpeg",path="../analysis/")

### ML-sequence
df = subset(read.csv("ml-sequence-N10000-L50.csv"), θ==0.1)
p <- df$P.bc>df$P.ab
q <- df$P.ac>df$P.ab
inconsistency_zone <- p | q
plot3 = ggplot(df, aes(x = ρ_a, y = τ_abc-τ_ab, color = inconsistency_zone)) +
    geom_point() + 
    labs(title="Approximate inconsistency zone under ML-sequence for θ=0.1",x="ρ_a",y="f=τ_abc-τ_ab")
ggsave("plot4-mls-th0.1.jpeg",path="../analysis/")

df = subset(read.csv("ml-sequence-N10000-L50.csv"), θ==0.01)
p <- df$P.bc>df$P.ab
q <- df$P.ac>df$P.ab
inconsistency_zone <- p | q
plot3 = ggplot(df, aes(x = ρ_a, y = τ_abc-τ_ab, color = inconsistency_zone)) +
    geom_point() + 
    labs(title="Approximate inconsistency zone under ML-sequence for θ=0.01",x="ρ_a",y="f=τ_abc-τ_ab")
ggsave("plot4-mls-th0.01.jpeg",path="../analysis/")

### JC-sequence
df = subset(read.csv("jc-sequence-N10000-L50.csv"), θ==0.1)
p <- df$P.bc>df$P.ab
q <- df$P.ac>df$P.ab
inconsistency_zone <- p | q 
plot3 = ggplot(df, aes(x = ρ_a, y = τ_abc-τ_ab, color = inconsistency_zone)) +
    geom_point() + 
    labs(title="Approximate inconsistency zone under JC-sequence for θ=0.1",x="ρ_a",y="f=τ_abc-τ_ab")
ggsave("plot4-jcs-th0.1.jpeg",path="../analysis/")

df = subset(read.csv("jc-sequence-N10000-L50.csv"), θ==0.01)
p <- df$P.bc>df$P.ab
q <- df$P.ac>df$P.ab
inconsistency_zone <- p | q
plot3 = ggplot(df, aes(x = ρ_a, y = τ_abc-τ_ab, color = inconsistency_zone)) +
    geom_point() + 
    labs(title="Approximate inconsistency zone under JC-sequence for θ=0.01",x="ρ_a",y="f=τ_abc-τ_ab")
ggsave("plot4-jcs-th0.01.jpeg",path="../analysis/")



##______________________________________________________________________________
##
## Part 5 -- Discrete illustrationss of inconsistency zone for inference of rooted
##           triplet topology.
##______________________________________________________________________________

## parameters same as for plot 3,4, but with L=500 rather than L=50
## made with datafiles in folder 2020-05-14

# JC expected color test
df = subset(read.csv("jc-expected-N10000-L500.csv"), θ==0.1)
p <- df$P.bc>df$P.ab
q <- df$P.ac>df$P.ab
inconsistency_zone <- p | q # the vertical bar | is a vectorized logical 'or'
plot3 = ggplot(df, aes(x = ρ_a, y = τ_abc-τ_ab, color = inconsistency_zone, size = P.bc - P.ab)) +
    geom_point() + 
    labs(title="Approximate inconsistency zone under JC-expected for θ=0.1,L=500",x="ρ_a",y="τ_AB")
ggsave("color-test-plot5-jce-th0.1.jpeg",path="../analysis/")

# JC-expected
df = subset(read.csv("jc-expected-N10000-L500.csv"), θ==0.1)
p <- df$P.bc>df$P.ab
q <- df$P.ac>df$P.ab
inconsistency_zone <- p | q # the vertical bar | is a vectorized logical 'or'
plot3 = ggplot(df, aes(x = ρ_a, y = τ_abc-τ_ab, color = inconsistency_zone)) +
    geom_point() + 
    labs(title="Approximate inconsistency zone under JC-expected for θ=0.1,L=500",x="ρ_a",y="f=τ_abc-τ_ab")
ggsave("plot5-jce-th0.1.jpeg",path="../analysis/")

df = subset(read.csv("jc-expected-N10000-L500.csv"), θ==0.01)
p <- df$P.bc>df$P.ab
q <- df$P.ac>df$P.ab
inconsistency_zone <- p | q
plot3 = ggplot(df, aes(x = ρ_a, y = τ_abc-τ_ab, color = inconsistency_zone)) +
    geom_point() + 
    labs(title="Approximate inconsistency zone under JC-expected for θ=0.01,L=500",x="ρ_a",y="f=τ_abc-τ_ab")
ggsave("plot5-jce-th0.01.jpeg",path="../analysis/")

### ML-sequence
df = subset(read.csv("ml-sequence-N10000-L500.csv"), θ==0.1)
p <- df$P.bc>df$P.ab
q <- df$P.ac>df$P.ab
inconsistency_zone <- p | q
plot3 = ggplot(df, aes(x = ρ_a, y = τ_abc-τ_ab, color = inconsistency_zone)) +
    geom_point() + 
    labs(title="Approximate inconsistency zone under ML-sequence for θ=0.1,L=500",x="ρ_a",y="f=τ_abc-τ_ab")
ggsave("plot5-mls-th0.1.jpeg",path="../analysis/")

df = subset(read.csv("ml-sequence-N10000-L500.csv"), θ==0.01)
p <- df$P.bc>df$P.ab
q <- df$P.ac>df$P.ab
inconsistency_zone <- p | q
plot3 = ggplot(df, aes(x = ρ_a, y = τ_abc-τ_ab, color = inconsistency_zone)) +
    geom_point() + 
    labs(title="Approximate inconsistency zone under ML-sequence for θ=0.01,L=500",x="ρ_a",y="f=τ_abc-τ_ab")
ggsave("plot5-mls-th0.01.jpeg",path="../analysis/")

### JC-sequence
df = subset(read.csv("jc-sequence-N10000-L500.csv"), θ==0.1)
p <- df$P.bc>df$P.ab
q <- df$P.ac>df$P.ab
inconsistency_zone <- p | q 
plot3 = ggplot(df, aes(x = ρ_a, y = τ_abc-τ_ab, color = inconsistency_zone)) +
    geom_point() + 
    labs(title="Approximate inconsistency zone under JC-sequence for θ=0.1,L=500",x="ρ_a",y="f=τ_abc-τ_ab")
ggsave("TEST-plot5-jcs-th0.1.jpeg",path="../analysis/")

df = subset(read.csv("jc-sequence-N10000-L500.csv"), θ==0.01)
p <- df$P.bc>df$P.ab
q <- df$P.ac>df$P.ab
inconsistency_zone <- p | q
plot3 = ggplot(df, aes(x = ρ_a, y = τ_abc-τ_ab, color = inconsistency_zone)) +
    geom_point() + 
    labs(title="Approximate inconsistency zone under JC-sequence for θ=0.01,L=500",x="ρ_a",y="f=τ_abc-τ_ab")
ggsave("TEST-plot5-jcs-th0.01.jpeg",path="../analysis/")




##______________________________________________________________________________
##
## Part 6 -- Effect of recombination when all other variables are fixed and
##          internal branch length is small and τ_ab varies
##______________________________________________________________________________

## Simulation parameters:
## (defparameter *τ_ab-values* '(.1 .5 1 2))
## (defparameter *f-values* '(.01))
## (defparameter *ρ_a-values* '(0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20))
## (defparameter *ρ_b-values* '(0))
## (defparameter *ρ_c-values* '(0))
## (defparameter *ρ_ab-values* '(0))
## (defparameter *ρ_abc-values* '(0))
## (defparameter *θ-values* '(.1)) 
## (defparameter *N* 15000)
## (defparameter *L* 500)
## (defparameter *τ_max* 999999)

## jc-expected
df = read.csv("jc-expected-N15000-L500.csv")
plot1 = ggplot(df, aes(x=ρ_a, y= P.bc - P.ab)) +
    geom_point() + 
    labs(title="Recombination effect on most likely inferred gene topology under JC-expected", x="recombination rate in population A",y="Pr[A(BC)]-Pr[(AB)C]")
ggsave("plot6-jce.jpeg",path="../analysis/")

## jc-sequence
df = read.csv("jc-sequence-N15000-L500.csv")
plot1jce = ggplot(df, aes(x=ρ_a, y= P.bc - P.ab)) +
    geom_point() + 
    labs(title="Recombination effect on most likely inferred gene topology under JC-sequence",x="recombination rate in population A",y="Pr[A(BC)]-Pr[(AB)C]")
ggsave("plot6-jcs.jpeg",path="../analysis/")

## ml-sequence
df = read.csv("ml-sequence-N15000-L500.csv")
plot1jce = ggplot(df, aes(x=ρ_a, y= P.bc - P.ab)) +
    geom_point() + 
    labs(title="Recombination effect on most likely inferred gene topology under ML-sequence",x="recombination rate in population A",y="Pr[A(BC)]-Pr[(AB)C]")
ggsave("plot6-mls.jpeg",path="../analysis/")




##______________________________________________________________________________
##
## Part 9 -- Recomb in A vs Recomb in B
##          
##______________________________________________________________________________

## Simulation parameters:
## Part 9. 
## τ_ab ∈ (1) 
## f ∈ (0.01) 
## τ_max=999999 
## θ ∈ (0.01) 
## N=10000 
## L=500 
## ρ_a ∈ (0 2 4 6 8 10 12 14 16 18 20) 
## ρ_b ∈ (0 2 4 6 8 10 12 14 16 18 20) 
## ρ_c ∈ (0) 
## ρ_ab ∈ (0) 
## ρ_abc ∈ (0)


## Part 9.1 -- Red/blue/black dot plot
df = read.csv("part9/part9-sim1-jc-sequence-N10000-L500.csv",header=TRUE)
x = df$ρ_a # x = recomb rate in population A
y = df$ρ_b # y = recomb rate in population B

p = ggplot(df, aes(x=ρ_a, y=ρ_b, size=2, color = factor(1*(P.bc>P.ab & P.bc>P.ac) - 1*(P.ac>P.ab & P.ac>P.bc) ))) + # categoral variable indicating which topology wins, with 1,2,and 3 corresponding to AC|B, AB|C, and BC|A respectively)
    geom_point() +
    scale_colour_manual(labels = c("AC|B", "AB|C", "AC|B"), values = c("black","#00BFC4FF","#F9766EFF")) +
    labs(title="Inference with Recombination in A and B",
         color="Most inferred topology",
         x=expression('Recombination rate in population A (ρ'['A'] * ')'),
         y=expression('Recombination rate in population B (ρ'['B'] * ')')) +
    scale_size(guide=FALSE) +
    theme(text=element_text(size=15)) # make text bigger
p
ggsave("part9-plot1.jpeg",path="../analysis/")



## Part 9.2 -- Red/blue/black dot plot
df = read.csv("part9/part9-sim2-jc-sequence-N10000-L500.csv",header=TRUE)
x = df$ρ_a # x = recomb rate in population A
y = df$ρ_b # y = recomb rate in population B

p = ggplot(df, aes(x=ρ_a, y=ρ_b, size=2, color = factor(1*(P.bc>P.ab & P.bc>P.ac) - 1*(P.ac>P.ab & P.ac>P.bc) ))) + # categoral variable indicating which topology wins, with 1,2,and 3 corresponding to AC|B, AB|C, and BC|A respectively)
    geom_point() +
    scale_colour_manual(labels = c("AC|B", "AB|C", "AC|B"), values = c("black","#00BFC4FF","#F9766EFF")) +
    labs(title="Inference with Recombination in A and B",
         color="Most inferred topology",
         x=expression('Recombination rate in population A (ρ'['A'] * ')'),
         y=expression('Recombination rate in population B (ρ'['B'] * ')')) +
    scale_size(guide=FALSE) +
    theme(text=element_text(size=15)) # make text bigger
p
ggsave("part9-plot2.jpeg",path="../analysis/")






## Part 9.3 -- 3D Surface plot
library(akima)
library(rgl)

rgl.open()
open3d()

## Interpolate
df = read.csv("/part9/part9-sim2-jc-sequence-N10000-L500.csv",header=TRUE)
x = df$ρ_a # x = recomb rate in population A
y = df$ρ_b # y = recomb rate in population B
z=df$P.ab


n_interpolation <- 100
spline_interpolated <- interp(x, y, z,
                              xo=seq(min(x), max(x), length = n_interpolation),
                              yo=seq(min(y), max(y), length = n_interpolation),
                              linear = TRUE, extrap = FALSE)
x.si <- spline_interpolated$x
y.si <- spline_interpolated$y
z.si <- spline_interpolated$z

## Number of colors + color vector
zcol  = cut(z.si, -1:1)
color = c("#F9766EFF","#00BFC4FF")


## Plot it
persp3d(x.si, y.si, z.si, col = color[zcol], ticktype="detailed", shade = 0.3, xlab = "", ylab = "", zlab = "")

decorate3d(main = "", xlab=expression(bold('ρ' ['A'])), ylab=expression(bold('ρ' ['B'])), zlab=expression(bold('z')))
rglwidget()






## Part 10. Redoing JC-sequence plots from part 5


### Plot 10.1. Cross-section of anomaly zone (red/blue dots)
### When this is repeated with θ=0.01 (as well as JCS) the result is similar.
df = subset(read.csv("part5/jc-sequence-N10000-L500.csv", header=TRUE), θ==0.1)
p <- df$P.bc>df$P.ab
q <- df$P.ac>df$P.ab
inconsistency_zone <- p | q 
plot = ggplot(df, aes(x = ρ_a, y = τ_abc-τ_ab, color = !inconsistency_zone, size =2)) +
    geom_point() + 
    labs(title="Inconsistency Zone of R*", #note θ=0.1,L=500
         x = expression('Recombination rate in population A (ρ' ['A'] * ')'),
         y= expression('Internal branch length (τ' ['AB'] * ')'),
         color="Correct inference?") +
    scale_size(guide=FALSE) + #removes the 'size' legend'
    guides(color = guide_legend(reverse = TRUE)) +
    theme(text=element_text(size=15)) # make text bigger
#print(plot)

plot

ggsave("part10-plot1.jpeg",path="../analysis/")






### Plot 10.2. 3D anomaly zone plot
### Got color code from https://stackoverflow.com/questions/17258787/formatting-of-persp3d-plot
### Adapted code from: https://stackoverflow.com/questions/36413652/3d-surface-interpolation?rq=1
library(akima)
library(rgl)

rgl.open()
open3d()
df = subset(read.csv("part5/jc-sequence-N10000-L500.csv", header=TRUE), θ==0.1)

x = df$ρ_a # x = recomb rate in population A
y = df$τ_abc - df$τ_ab # y = internal branch length
z = df$P.ab - df$P.bc # z = P[ab|c]-P[bc|a]


# Interpolate
n_interpolation <- 25
spline_interpolated <- interp(x, y, z,
                              xo=seq(min(x), max(x), length = n_interpolation),
                              yo=seq(min(y), max(y), length = n_interpolation),
                              linear = TRUE, extrap = FALSE)
x.si <- spline_interpolated$x
y.si <- spline_interpolated$y
z.si <- spline_interpolated$z


# Number of colors + color vector
#nbcol = 100
#color = rainbow(nbcol, start = 0/6, end = 4/6)
#color
zcol  = cut(z.si, -1:1)
color = c("#F9766EFF","#00BFC4FF")
color


## Plot it
persp3d(x.si, y.si, z.si, col = color[zcol], ticktype="detailed", shade = 0.3, xlab = "", ylab = "", zlab = "")
decorate3d(main = "", xlab=expression(bold('ρ' ['A'])), ylab=expression(bold('τ' ['AB'])), zlab=expression(bold('z')))
rglwidget()






## Part 11. All the same recombination rates

## (defparameter *τ_ab-values* '(1))
## (defparameter *f-values* '(.01))
## ρ_a, ρ_b ρ_c ρ_ab were all equal for this part, and took values (1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20)
## (defparameter *ρ_abc-values* '(0))
## (defparameter *θ-values* '(.01))
## (defparameter *N* 10000)
## (defparameter *L* 500)



getwd()
df=read.csv("part11/part11-sim2-jc-sequence-N100000-L500.csv")


p = ggplot(df, aes(x=ρ_a, y=P.ab-pmax(P.ac,P.bc), size=5, color = TRUE)) + 
    geom_point() +
    scale_colour_manual(values = c("#00BFC4FF")) +
    labs(title="Equal Recombination Rates",
         x=expression('Recombination rate ρ'),
         y=expression('P[AB|C]-max(P[AC|B],P[BC|A])')) +
    scale_size(guide=FALSE) + # remove scale legend
    guides(color=FALSE) + # remove color legend
    theme(text=element_text(size=20)) # make text bigger

p

ggsave("part11-plot.jpeg",path="../analysis/")





