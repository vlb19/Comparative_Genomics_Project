#!/usr/bin/env Rscript

# Author: Victoria Blanchard vlb19@ic.ac.uk
# Script: WrightFisherTutorial.R
# Desc: Demonstrate	how	the Wright-Fisher model works
# Input: 
# Output: 
# Arguments: 0
# Date: March 2020
# Requirements: Ape package (sudo apt-get install -y -r-cran-ape)

##############################################
### Prepare R environment ###
##############################################

#Clear global environment
rm(list=ls()) 

##############################################
### Source Wright-Fisher model script ###
##############################################

source("simulateWF.R")

### This script contains a set of functions for simulating the Wright-Fisher model, both forwards and backwards in time

##############################################
## 1 - Thinking forwards in time: 2 alleles
##############################################

### Run a Wright-Fisher model beginning a population with two allels
# The population will have size 2N = 5 (so N = 5 diploids) and we'll 
# run the simulation for 15 generations

WF_twoalleles(5, 15)

### Save the vector of blue allele counts and trace the frequency of 
# the blue allele over time 

bluecounts <- WF_twoalleles(5, 15)
bluefreq <- bluecounts / (2*5)
plot(bluefreq, ylim=c(0,1), type = "b", col = "blue", pch=19, xlab = "generations", ylab = "Blue frequency")

### Do alleles tend to "fix" faster when N is large or when N is small?

"They tend to fix faster when N is small"

##############################################
## 2 - Thinking forwards in time: many alleles
##############################################

### We can aslo run a Wright-Fisher model with more than two alleles
# In this simulation, each individual contains two distinct alleles 
# which are different from all other alleles in the population 

WF_manyalleles(5, 15)

# What happens to the allelic diversity as time goes forward?

"It drastically reduces over time"

# Are there more or less heterozygotes at the end of the simulation than 
# at the beginning?

"Many of the original alleles become extinct almost immediately. There 
are mainly homozygotes at the end generally"

# What happens to allelic diversity over time, when N = 3 and N = 10?

"Allelic diversity diminishes almost immediately for low values of N, 
whereas a much larger proportion of the alleles do not go extinct in 
large values of N"

##############################################
## 3 - Thinking backwards in time
##############################################

### So far these simulations have been running forwards in time. We will 
# now trace the lineages of particular individuals from the present 
# backwards in time to see how they coalesce (find a common ancestor) in 
# the past

### We will trace the genealogy of 3 lineages in a population of size 
# N = 10 (2N = 20) over 20 generations:

track_lineages(N.vec = rep(10,20), n.iter = 1, num.tracked = 3)

### Repeat this simulation 10 times. For each simulation record the time 
# between the present and the first coalescent event, and the time between 
# the first coalescent event and the second coalescent event (i.e. the
# most recent common ancesotr of all 3 lineages). 

"1) 6 generations, 11 generations
 2) 2 generations, 16 generations
 3) 12 generations, 6 generations
 4) 7 generations, 8 generations 
 5) 14 generations, 3 generations 
 6) 7 generations, 6 generations 
 7) 1 generation, 17 generations 
 8) 10 generations, 5 generations 
 9) 2 generations, 12 generations 
10) 15 generations, 1 generation "

# Which of the two times tends to be larger? Why?

"10 repeats did not reveal this but overall it seems like the time from 
the first coalescence event to the second seems to be longer. This is 
because as we go backwards in time, we will have less and less sequences
so the time to the next coalescent event will be larger."

# Check what happens to the coalescence rate when N = 7 and when N = 20
# Do lineages coalesce faster or slower with larger population size?

track_lineages(N.vec = rep(7,14), n.iter = 1, num.tracked = 3)
track_lineages(N.vec = rep(20,40), n.iter = 1, num.tracked = 3)

"Lineages coalesce faster with increasing values of N. "