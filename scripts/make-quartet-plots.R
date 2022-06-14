##______________________________________________________________________________
##
## Setup -- Install and load necessary R packages.
##______________________________________________________________________________

## You need to run the code in this section before doing anything else. 

packages <- c("ggplot2","gridExtra")

installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages]) # Install packages
}

invisible(lapply(packages, library, character.only = TRUE)) # Packages loading

setwd('../data/') # Set working directory to 'data/' 





##______________________________________________________________________________
##
## Quartet Experiment 1
##______________________________________________________________________________

df = read.csv("quartet-simulation-1.csv", header=TRUE)

zone <- as.factor(1*(df$AD.BC > df$AB.CD & df$AD.BC > df$AC.BD) + 2*(df$AC.BD > df$AB.CD & df$AC.BD > df$AD.BC))
                                        # categorical variable indicating which
                                        # quartet topology wins, with 0,1,and 2
                                        # corresponding to AB|CD, AD|BC, and
                                        # AC|BD respectively)

top_names <- c(expression(hat('q')*'=AB|CD'), expression(hat('q')*'=AD|BC'),  expression(hat('q')*'=AC|BD'))

p = ggplot(df, aes(x=ρ, y=f.one, size=2, group=zone))+
    geom_point(aes(shape=zone, color=zone, size=2))+
    scale_shape_manual(labels=top_names, values=c(16,17,15))+
    scale_colour_manual(labels=top_names, values = c("#00BFC4FF","#F9766EFF","black"))+
    labs(color="Uniquely favored\nunrooted quartet",
         shape="Uniquely favored\nunrooted quartet",
         x=expression('Recombination rate in populations A and B'),
         y=expression('Internal branch length')) +     
    scale_size(guide=FALSE) + # omit size on legend
    theme(text=element_text(size=20), legend.text=element_text(size=20))+ # make text bigger
    scale_y_continuous(limits=c(0.005, 0.065), expand=c(0,0))

p +  guides(color = guide_legend(override.aes = list(size=4))) -> p # make legend text bigger

p

ggsave("quartet-simulation-1.jpeg",path="../analysis/")




##______________________________________________________________________________
##
## Quartet Experiment 2
##______________________________________________________________________________


df=read.csv("quartet-simulation-2.csv", header=TRUE)


p = ggplot(df, aes(x=ρ, y=(AB.CD-pmax(AC.BD,AD.BC))/m, size=5, color = TRUE)) + 
    geom_point() +
    scale_colour_manual(values = c("#00BFC4FF")) +
    labs(#title="Equal Recombination Rates",
        x=expression('ρ'),
        y=expression(widehat('p')['AB|CD']*' - max('*widehat(p)['AC|BD']*' , '*widehat(p)['AD|BC']*')')) +
    scale_size(guide=FALSE) + # remove scale legend
    guides(color=FALSE) + # remove color legend
    theme(text=element_text(size=20)) + # make text bigger
    ylim(0, 0.02)

ggsave("quartet-simulation-2-plot.jpeg",path="../analysis/")



##______________________________________________________________________________
##
## Quartet Experiment 3
##______________________________________________________________________________

df = read.csv("quartet-simulation-3.csv", header=TRUE)

zone <- as.factor(1*(df$AD.BC > df$AB.CD & df$AD.BC > df$AC.BD) + 2*(df$AC.BD > df$AB.CD & df$AC.BD > df$AD.BC))
                                        # categorical variable indicating which
                                        # quartet topology wins, with 0,1,and 2
                                        # corresponding to AB|CD, AD|BC, and
                                        # AC|BD respectively)

top_names <- c(expression(hat('q')*'=AB|CD'), expression(hat('q')*'=AD|BC'),  expression(hat('q')*'=AC|BD'))

p = ggplot(df, aes(x=ρ, y=f, size=2, group=zone))+
    geom_point(aes(shape=zone, color=zone, size=2))+
    scale_shape_manual(labels=top_names, values=c(16,17,15))+
    scale_colour_manual(labels=top_names, values = c("#00BFC4FF","#F9766EFF","black"))+
    labs(color="Uniquely favored\nunrooted quartet",
         shape="Uniquely favored\nunrooted quartet",
         x=expression('ρ'),
         y=expression('Internal branch length')) +     
    scale_size(guide=FALSE) + # omit size on legend
    theme(text=element_text(size=20), legend.text=element_text(size=20))+ # make text bigger
    scale_y_continuous(limits=c(0.005, 0.105), expand=c(0,0))

p +  guides(color = guide_legend(override.aes = list(size=4))) -> p # make legend text bigger

ggsave("quartet-simulation-3.jpeg",path="../analysis/")


##______________________________________________________________________________
##
## Quartet Experiment 4
##______________________________________________________________________________

df=read.csv("quartet-simulation-4.csv", header=TRUE)

p = ggplot(df, aes(x=ρ, y=(AB.CD-pmax(AC.BD,AD.BC))/m, size=5, color = TRUE)) + 
    geom_point() +
    scale_colour_manual(values = c("#00BFC4FF")) +
    labs(#title="Equal Recombination Rates",
        x=expression('ρ'),
        y=expression(widehat('p')['AB|CD']*' - max('*widehat(p)['AC|BD']*' , '*widehat(p)['AD|BC']*')')) +
    scale_size(guide=FALSE) + # remove scale legend
    guides(color=FALSE) + # remove color legend
    theme(text=element_text(size=20)) + # make text bigger
    ylim(0.01, 0.025)

ggsave("quartet-simulation-4-plot.jpeg",path="../analysis/")


##______________________________________________________________________________
##
## Quartet Experiment 5
##______________________________________________________________________________

df = read.csv("quartet-simulation-5.csv", header=TRUE)

zone <- as.factor(1*(df$AD.BC > df$AB.CD & df$AD.BC > df$AC.BD) + 2*(df$AC.BD > df$AB.CD & df$AC.BD > df$AD.BC))
                                        # categorical variable indicating which
                                        # quartet topology wins, with 0,1,and 2
                                        # corresponding to AB|CD, AD|BC, and
                                        # AC|BD respectively)

top_names <- c(expression(hat('q')*'=AB|CD'), expression(hat('q')*'=AD|BC'),  expression(hat('q')*'=AC|BD'))

p = ggplot(df, aes(x=ρ, y=f, size=2, group=zone))+
    geom_point(aes(shape=zone, color=zone, size=2))+
    scale_shape_manual(labels=top_names, values=c(16,17,15))+
    scale_colour_manual(labels=top_names, values = c("#00BFC4FF","#F9766EFF","black"))+
    labs(color="Uniquely favored\nunrooted quartet",
         shape="Uniquely favored\nunrooted quartet",
         x=expression('ρ'),
         y=expression('Internal branch length')) +     
    scale_size(guide=FALSE) + # omit size on legend
    theme(text=element_text(size=20), legend.text=element_text(size=20))+ # make text bigger
    scale_y_continuous(limits=c(0.005, 0.105), expand=c(0,0))

p +  guides(color = guide_legend(override.aes = list(size=4))) -> p # make legend text bigger

ggsave("quartet-simulation-5.jpeg",path="../analysis/")
