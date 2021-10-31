## 2021-10-28 R Intro
## Source: http://leg.ufpr.br/~silvia/R2001/contents.html
## This document contains my random experiments with R for making plots.


## Store values
x <- 3
y <- x
x <- 2

## Vectors
v <- 1:10
v
y <- c( v, 2, 3,4)


z <- c(1:3, 1:4)
z


z= replicate(420,0)
z

## Define a (finite) sequence with a given start, end, and step size
seq(1, 10, .5)
seq(10, 1, -1)
seq(10, 1, -2)
seq(3, 87)

## Vector Indexes
v <- 1:10
v+2 -> v
v[1]

v[1:3]
v[100000]
v[-1] # nice way to remove an element from a vector
v[-3] # works for any index number :)
v
v==4
v == 4
v=3
v
w=rev(v) # rev() reverses the vector
w

## Managing variables 
ls()
rm(x)
ls()

## Useful functions
summary(v)


## Loading csv
x = read.csv("foo.csv")
x[1]
x[2,3]

## Plotting
x <- 1:20
y <- x**3
plot(x,y)
lines(x,y)
points(rev(x),y,pch=3)

points(x,y,pch="$")

lines(x,y,lwd=4) # thicc line
lines(rev(x),y,lty=2) # dash line
lines(rev(x),y,lwd=5,lty=2) # thicc dash line


## Print plot to a file
dev.print(file="gme-stock-price.ps")
plot(x,y,xlab="Monke hold banana long tim",ylab="Rocket Up", pch="$")
dev.off()





## For loops
for (year in 2010:2015){
  print(paste("The year is", year))
}



## Automate plot creation
v=1:10

for (item in v)
{
plot(x,y,xlab="Monke hold banana long tim",ylab="Rocket Up", pch="$")
dev.print(file=sprintf("gme-stonk-%s.ps", item)) 
dev.off()
}

## In the above, sprintf can accomodate more variables. Use the notation sprintf(“%s_%s.csv”,var1, var2)
a
b=2
sprintf("%s...%s",a,b)
sprintf("%s_%s",a,b)

## Exit R
q()

## To run R from shell
Rscript -e '1+1'


## First steps: Connecting the pipeline (2021-03-09)
##
## First load simulator.fasl and execute-consensus-jc.fasl. Then from the
## scripts directory, run the command:
##
## (mapcar #'try-various-recombination-rates '(160000 80000 40000 20000 10000) '(2 4 8 16 32))
##
## This generates fives files, located in the data/ directory:
##
##   consensus-jc-N10000-L32-F0.1.csv
##   consensus-jc-N20000-L16-F0.1.csv
##   consensus-jc-N40000-L8-F0.1.csv
##   consensus-jc-N80000-L4-F0.1.csv
##   consensus-jc-N160000-L2-F0.1.csv
##   
## Each of such file has 32 rows, one row for every possible assignment of
## recombination rate 0 or 1 to the five populations in the species tree (for
## example (1 0 0 0 0) or (0 1 1 0 1) or whatever). In particular, each row is
## generated by simulating N=10000,20000,40000,...,160000 ancestral
## recombination graphs for loci of size L=32,16,...,4,2 respectively. Loci of
## shorter length require more samples to get an accurate estimate of the
## inference probabilities. I haven't checked if this scaling is sufficient, but
## I'm short on time right now (2021-03-09) and I'm still figuring the pipeline out.
##
## Next, load up the R environment. To load these csv files as dataframes into
## R, run the commands

x1 = read.csv("consensus-jc-N160000-L2-F0.1.csv",header=TRUE)
x2 = read.csv("consensus-jc-N80000-L4-F0.1.csv)",header=TRUE)
x3 = read.csv("consensus-jc-N40000-L8-F0.1.csv",header=TRUE)
x4 = read.csv("consensus-jc-N20000-L16-F0.1.csv",header=TRUE)
x5 = read.csv("consensus-jc-N10000-L32-F0.1.csv",header=TRUE)

## (The header=TRUE is not required, since this is the default.)
##
## Now we want to analyze the data and make some pictures. Useful source:
## https://www.stat.berkeley.edu/~s133/R-4a.html
##
## The estimated probability of inferring the true species topology:

summary(x1[1])
summary(x2[1])
summary(x3[1])
summary(x4[1])
summary(x5[1])

## and we can do things like this, which doesn't look like much

plot(1:32,x5$P.ab)

## and this, which seems to suggest we better increase our samples

boxplot(x1$P.ab, x2$P.ab, x3$P.ab, x4$P.ab, x5$P.ab)




## 2021-03-10
## in terminal, navigate to data/data-correctly-named/

cat consensus-jc-N20000-L32-F0.01-TH-0.1.csv | column -s, -t | 

## just look at cases where the true topology is less likely than P-ac
    cat consensus-jc-N20000-L32-F0.01-TH-0.1.csv | column -s, -t | awk '$1 < $2' | wc -l

## just look at cases where the true topology is less likely than P-bc
cat consensus-jc-N20000-L32-F0.01-TH-0.1.csv | column -s, -t | awk '$1 < $3' | wc -l

## For f=0.01, we have P-ab < min(P-ac, P-bc) in 114 out of 626 regimes
cat consensus-jc-N20000-L32-F0.01-TH-0.1.csv | awk -F',' '$1 < $3 && $1 < $2' | wc -l

## For f=0.1, apparent inconsistency occurs in 0 out of 626 regimes.

## Let's make a csv file with the candidates for inconsistency regimes
cat consensus-jc-N20000-L32-F0.01-TH-0.1.csv | awk -F',' '$1 < $3 && $1 < $2' > jc-inconsistency-candidates.csv
## To analyze this file, we can compute the averages of the recombination rates

for column in {7..11}
do
    awk -v column=$column -F',' '{sum+=$column} END{print sum/114}' jc-inconsistency-candidates.csv
done

4.65789
4.72807
0.675439
3.14035
0

# which gives the values of ρ_a,ρ_b,ρ_c,ρ_ab,ρ_abc. Note that ρ_abc=0 for all loci, so that last number doesn't tell us anything. Also, ρ_a >0 for all these 115 cases.


## Back in R, we can do things like: 
ic=read.csv("jc-inconsistency-candidates.csv")
summary(ic$P.ab)
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
0.3028  0.3148  0.3194  0.3191  0.3239  0.3309

## Look at the differences
> summary(ic$P.ac - ic$P.ab)
    Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
0.000400 0.008013 0.018650 0.021664 0.033063 0.062500 
> summary(ic$P.bc - ic$P.ab)
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
0.00090 0.00765 0.01773 0.02105 0.03204 0.05980

## Here are two interesting plots!  
plot(ic$ρ_a,
     ic$P.bc - ic$P.ab,
     xlab="recombination-rate_A",ylab="P.bc-P.ab", pch=1)
dev.print(file="../../analysis/plot-1-2021-03-12.ps") 
dev.off()

plot(ic$ρ_a - ic$ρ_c, ic$P.bc - ic$P.ab,xlab="recombination-rate_A - recombination-rate_C",ylab="P.bc-P.ab", pch=1)
dev.print(file="../../analysis/plot-2-2021-03-12.ps") 
dev.off()

## evidence this effect holds for ML as wll, but the effect is weaker:
(time (ml-estimate-topology-probabilities 1 1.01 999999 10 0 0 0 0 .1 1000000 32))
Evaluation took:
  123.135 seconds of real time
  123.379030 seconds of total run time (123.098681 user, 0.280349 system)
  357,588,545,367 processor cycles
  45,989,591,168 bytes consed
  
"0.33116034,0.32483435,0.34400535,1,1.01,999999,10,0,0,0,0,0.1,1000000,32"

## similar plots (need to fix labels)
plot(ic$ρ_b - ic$ρ_c, ic$P.bc - ic$P.ac,xlab="recombination-rate_A - recombination-rate_C",ylab="P.bc-P.ab", pch=1)
dev.print(file="../../analysis/plot-2-2021-03-12.ps") 
dev.off()

plot(x2$ρ_a,
     x2$P.bc - x2$P.ab,
     xlab="recombination-rate_A",ylab="P.bc-P.ab", pch=1)
dev.print(file="../../analysis/plot-0-2021-03-12.ps") 

## Check out these parameter values
(time (ml-estimate-topology-probabilities 1 1.01 999999 10 0 0 0 0 .02 1000 2000))

Evaluation took:
  141.393 seconds of real time
  141.397203 seconds of total run time (141.393376 user, 0.003827 system)
  [ Run times consist of 0.066 seconds GC time, and 141.332 seconds non-GC time. ]
  410,612,929,727 processor cycles
  1,701,198,288 bytes consed
  
"0.312,0.2965,0.3915,1,1.01,999999,10,0,0,0,0,0.02,1000,2000"

## and again
(time (estimate-topology-probabilities 1 1.01 999999 10 0 0 0 0 .02 1000 2000))
Evaluation took:
  120.048 seconds of real time
  120.051731 seconds of total run time (120.015893 user, 0.035838 system)
  [ Run times consist of 0.068 seconds GC time, and 119.984 seconds non-GC time. ]
  348,623,048,955 processor cycles
  1,690,918,352 bytes consed
  
"0.313,0.296,0.391,1,1.01,999999,10,0,0,0,0,0.02,1000,2000"


## again
(time (estimate-topology-probabilities 1 1.01 999999 10 0 0 0 0 .002 1000 2000))
Evaluation took:
  119.280 seconds of real time
  119.274739 seconds of total run time (119.174731 user, 0.100008 system)
  [ Run times consist of 0.070 seconds GC time, and 119.205 seconds non-GC time. ]
  346,392,421,085 processor cycles
  1,687,874,176 bytes consed
  
"0.319,0.305,0.376,1,1.01,999999,10,0,0,0,0,0.002,1000,2000"

# and again
(time (estimate-topology-probabilities 1 1.01 999999 10 0 0 0 0 .02 1000 2000))
Evaluation took:
  122.801 seconds of real time
  122.817544 seconds of total run time (122.813585 user, 0.003959 system)
  [ Run times consist of 0.070 seconds GC time, and 122.748 seconds non-GC time. ]
  356,621,412,702 processor cycles
  1,684,325,616 bytes consed
  
"0.309,0.297,0.394,1,1.01,999999,10,0,0,0,0,0.02,1000,2000"




## Part 10. Redoing JC-sequence plots from part 5

packages <- c("ggplot2")
invisible(lapply(packages, library, character.only = TRUE)) # Packages loading
setwd('../scripts/')


### Plot 10.1. Cross-section of anomaly zone (red/blue dots)
### When this is repeated with θ=0.01 (as well as JCS) the result is similar.
df = subset(read.csv("../data/part5/jc-sequence-N10000-L500.csv", header=TRUE), θ==0.1)
p <- df$P.bc>df$P.ab
q <- df$P.ac>df$P.ab
inconsistency_zone <- p | q 
plot = ggplot(df, aes(x = ρ_a, y = τ_abc-τ_ab, color = !inconsistency_zone, size =2)) +
    geom_point() + 
    labs(title="Inconsistency Zone of R* with sequence distances", #note θ=0.1,L=500
         x = expression('Recombination rate in population A (ρ' ['A'] * ')'),
         y= expression('Internal branch length (τ' ['AB'] * ')'),
         color="Topology correctly inferred?") +
    scale_size(guide=FALSE) + #removes the 'size' legend'
    guides(color = guide_legend(reverse = TRUE))
#print(plot)
ggsave("part10-plot1.jpeg",path="../analysis/")






### Plot 10.2. 3D anomaly zone plot
### Got color code from https://stackoverflow.com/questions/17258787/formatting-of-persp3d-plot
### Adapted code from: https://stackoverflow.com/questions/36413652/3d-surface-interpolation?rq=1
rgl.open()
open3d()
df = subset(read.csv("../data/part5/jc-sequence-N10000-L500.csv", header=TRUE), θ==0.1)

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









## Part 9. Recomb in A vs Recomb in B
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

##> getwd()
##[1] "/home/mutalisk/MSCR-simulator/scripts"


## Part 9.1 -- Red/blue/black dot plot
packages <- c("ggplot2")
invisible(lapply(packages, library, character.only = TRUE)) # Packages loading
setwd('../data/part9/')

df = read.csv("part9-sim1-jc-sequence-N10000-L500.csv",header=TRUE)
x = df$ρ_a # x = recomb rate in population A
y = df$ρ_b # y = recomb rate in population B

p = ggplot(df, aes(x=ρ_a, y=ρ_b, size=2, color = factor(1*(P.bc>P.ab & P.bc>P.ac) - 1*(P.ac>P.ab & P.ac>P.bc) ))) + # categoral variable indicating which topology wins, with 1,2,and 3 corresponding to AC|B, AB|C, and BC|A respectively)
    geom_point() +
    scale_colour_manual(labels = c("AC|B", "AB|C", "AC|B"), values = c("black","#00BFC4FF","#F9766EFF")) +
    labs(title="Inference with Recombination in A and B",
         color="Most inferred topology",
         x=expression('Recombination rate in population A (ρ'['A'] * ')'),
         y=expression('Recombination rate in population B (ρ'['B'] * ')')) +
    scale_size(guide=FALSE)
p
ggsave("part9-plot1.jpeg",path="../../analysis/")




ggsave("part9-plot1.jpeg",path="../analysis/")


## Part 9.3 -- 3D Surface plot
# Interpolate
df = read.csv("part9-sim2-jc-sequence-N10000-L500.csv",header=TRUE)
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

x.si
y.si
z.si

# Number of colors + color vector
#nbcol = 100
#color = rainbow(nbcol, start = 0/6, end = 4/6)
#color
zcol  = cut(z.si, -1:1)
color = c("#F9766EFF","#00BFC4FF")
color


## Plot it
persp3d(x.si, y.si, z.si, col = color[zcol], ticktype="detailed", shade = 0.3, xlab = "", ylab = "", zlab = "")

decorate3d(main = "", xlab=expression(bold('ρ' ['A'])), ylab=expression(bold('ρ' ['B'])), zlab=expression(bold('z')))
rglwidget()





## Part 11. All the same recombination rates

## (defparameter *τ_ab-values* '(1)) ; age of most recent species divergence

## (defparameter *f-values* '(.01))
## 				 ; f is the internal branch length on the
## 				 ; species tree: f=τ_abc-τ_ab.

## ;; The following parameters define population-specific recombination rates.
## ;; There are 5 populations, corresponding to edges on the species tree (A, B, C,
## ;; AB) as well as the root of the species tree (ABC). The simulator will
## ;; consider all possible combinations of the listed values.
## (defparameter *ρ_a-values* '(2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20))
## (defparameter *ρ_b-values* '(0))
## (defparameter *ρ_c-values* '(0))
## (defparameter *ρ_ab-values* '(0))
## (defparameter *ρ_abc-values* '(0))  ; Note that choosing a non-zero values for
## 				    ; ρ_abc can significantly increase computing
## 				    ; time.

## (defparameter *θ-values* '(.01)) ; mutation rate

## (defparameter *N* 10000) ; sample size (number of sampled loci)
## (defparameter *L* 500) ; locus length (in base pairs)


## (defparameter *τ_max* 999999)


getwd()
df=read.csv("../data/part11/part11-sim2-jc-sequence-N100000-L500.csv")



df$P.ab > df$P.ac & df$P.ab > df$P.bc

x=df$ρ_a
x
y=df$P.ab-pmax(df$P.ac,df$P.bc)
plot(x,y)



p = ggplot(df, aes(x=ρ_a, y=df$P.ab-pmax(df$P.ac,df$P.bc), size=2, color = TRUE)) + 
    geom_point() +
    scale_colour_manual(values = c("#F9766EFF")) +
    labs(title="Equal recombination in populations A, B, C, AB",
         x=expression('Recombination rate ρ in populations A, B, C, and AB'),
         y=expression('P[AB|C]-max(P[AC|B],P[BC|A]')) +
    scale_size(guide=FALSE)

p


