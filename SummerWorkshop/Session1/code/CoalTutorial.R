#!/usr/bin/env Rscript

# Author: Victoria Blanchard vlb19@ic.ac.uk
# Script: CoalTutorial.R
# Desc: Demonstrate	how	a	coalescence tree	is	simulated
# Input: 0
# Output: 0
# Arguments: 0
# Date: January 2020
# Requirements: Ape package (sudo apt-get install -y -r-cran-ape)

##############################################
### Prepare R environment ###
##############################################

#Clear global environment
rm(list=ls()) 

"1.) Start	by	drawing	on a piece of paper a	small circle	for	each	of	the	five	gene	copies. 
They should be lined up on	an	invisible horizontal line	and you should leave enough	space	above the circles for	
drawing	a	tree	above them (which we will do shortly).	We will henceforth call these five circles nodes and label them 1,2,3,4,5"

### 2.) Make a list of node names 

nodes = c(1,2,3,4,5) # make the list and call it nodes
nodes # print the list

### 3.) Sample	which	two	nodes	will	coalesce	first	(going	back	in	time)	by	randomly	picking	two	of	the	
### nodes.

nodecount = length(nodes) # save the number of nodes in the variable nodecount
tocoalesce = sample(1:nodecount, size=2) # sample 2 different nodes in node list
nodes[tocoalesce[1]] # print the first node sampled
nodes[tocoalesce[2]] # print the second node sampled

### 4.) Sample	the	time	it	takes	before	these	two	nodes	coalesce	(measured	from	previous	
### coalescence	event in	units	of	2N)	by	sampling	from	an	exponential	distribution	with	rate	equal	
### to	nodecount*(nodecount-1)/2	where	nodecount	is	the	number	of	nodes	in	your node	list.	

coalescencerate = nodecount*(nodecount-1)/2 # calculate the coalescent rate
coalescencetime = rexp(1, rate=coalescencerate) # sample from exponential with that rate 
coalescencetime # print coalescence time 

" 5.) Now	draw	a	node	that	is	the	sampled amount	of	time	further	up	in	the	tree	than	the currently	
highest node	(so	if	the	currently	highest	node	is	drawn	at	height	T	then	draw	the	new	one	at	
height	T plus the	sampled	coalescence	time)	and	draw	a	branch	from	each	of	the	nodes	you	
sampled	in	step	3	to	this	new	node	indicating	that	these	two	nodes	coalesce at	this	time."

### 6.) Next,	make	an	updated	list	of	the	nodes	that	are	left	by	removing	the two	nodes	that	
### coalesced	and	instead	adding	the	newly	drawn	node	that represents	their	common	ancestor.	

nodes <- nodes[-tocoalesce] # remove the two nodes that coalesced
nodes <- c(nodes,2*5-length(nodes)-1) # add the new node
nodes # print the new list

### 7.) If	you	only	have	one	node	left	in	your	list	of	remaining	nodes	you	are	done.	If	not,	go	back	to	
### step	3.	
for (i in 1:100) {
  if (length(nodes)>1){
    # Sample 2 different nodes in node list
    nodecount = length(nodes) # save the number of nodes in the variable nodecount
    tocoalesce = sample(1:nodecount, size=2) # sample 2 different nodes in node list
    
    # Sample time taken before coalescence
    coalescencerate = nodecount*(nodecount-1)/2 # calculate the coalescent rate
    coalescencetime = rexp(1, rate=coalescencerate) # sample from exponential with that rate 
    
    # Update the list of nodes 
    nodes <- nodes[-tocoalesce] # remove the two nodes that coalesced
    nodes <- c(nodes,2*5-length(nodes)-1) # add the new node
  }
}

nodes


### Doing	this	by	hand	is	obviously a	bit	tedious.	So	based	on	the	R	code	snippets	you	already	
### got, we have built a function	that	allows	you	to	do	this	automatically	(it	even	makes	a	drawing	of	the	tree).

source("simulatecoalescencetrees.R")

### Simulate	and	draw 10 trees just like	you	just like you did	by hand by typing the code below

par (mfrow=c(2,5))
for (i in c(1:10)){
  print("New Tree")
  yourtree <-simtree(5) # simulate tree with 5 nodes
  ct<-read.tree(text=yourtree);plot(ct,cex=1.5);add.scale.bar(y=1.2,x=0.2,cex = 2,col = "red",lcol="red",lwd=3)
  print(" ")
}


### 1) Which	coalescence event takes	the	longest on	average (the	first coalescence event,	the	
### second,	…,	or	the	last)?	And	which	event	takes	the	shortest on	average?

"The last coalescence event seems to take the longest on average and the first and second events 
seem to take the shortest amount of time" 
  
### 2) Is	that	what	you	would	expect? Recall that the	mean	of	an	exponential	distribution	with rate	lambda	
### is	1/lambda	and	the	coalescence rate	when	there	are	n nodes	left	is	n(n-1)/2.	So	the	mean	is	2/(n(n-1)),	so
### for	instance	for	when	there	are	5	nodes	left	the	mean	coalescent	time	is	2/(5(5-1))=0.1

"Yes it is to be expected, as n decreases, the mean coalescent time increases so the time to the last coalescence event is the longest"

### 3) Which	coalescence event	time	seems to	vary	the	most?

"The second coalescence event"
  
### 4) Is	that	what	you	would	expect? Recall that, if we have a random variable that follows an exponential distribution
### with rate lambda, then its variance is equal 1/(lambda^2).

"Yes it is to be expected "

"5) Finally,	imagine	the	following	case:	a	researcher	has	estimated	the	structure	of	a	tree	for	
mtDNA	from	a	species	sampled	in	a	single	location.	She	obtains	a	tree	looking	as	follows:
  
  ![alt text](https://github.com/FerRacimo/CopenhagenTutorial/blob/master/Tree0.png)

Based	on	the	structure	of	the	tree,	i.e.	two	groups	of	related	individuals	separated	by	long	
branches	down	to the	root	of	the	tree,	she	concludes	that	there	must	be population	
subdivision	with	two	clearly	differentiated	groups.	 Based	on	what	you	have	learned	from	
the	simulations,	do	you	agree	with	this	conclusion?"

# Yes because the branches of the tree are not sequential. There are multiple coalescence events happening
# at the same time and the products of those events are coalescing 
  