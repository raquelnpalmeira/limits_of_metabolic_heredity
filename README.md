# Code used for "The limits of metabolic heredity"

The files in this repository contain the necessary scripts to get the data presented in each of the figures in the paper.



## Figure 2


**partFunctionNucCat.m** is a function containing the main system of ODEs that form the model. It takes a 7 element vector x (containing number of molecules of C2, energy, fatty acids, and amino acids 1 and 2, sugars and nucleotides respectively) as input and returns a 7 element vector dx (containing the rate of change in the amount of the 7 species). 

#

**intFunctionNucCat.m** contains a function that performs the Euler integration of the system of ODE's in partFunctionBasic. The function takes three input arguments: K_C2_N, the rate constant for the catalysis of carbon fixation by nucleotides, ndays, a numerical value of the number of days the simulation should run for; and div, an argument that if equal to 1 allows for cell division. This function returns the average number of protocell divisions per day after 10 days of simulation. 

# 

Running the script Fig2 results in Figure 2. The script runs the function intFunctionNucCat.m with a range of values of K_C2_N are plots the vector with the results. 

## Figure 3

**partFunctionNucCatAA.m** and intFunctionNucCatAA.m respectively contain the system of ODEs and the Euler integration function for the formulation of the protocell model in which nucleotides catalyse carbon fixation and amino acid production, partFunctionNucCatAA.m takes in two extra inputs: K_AA_N, the rate constant for the catalysis of amino acid production by nucleotides and cost, a logical value for whether the model considers a cost to catalysis. Functions of the same name but ending in E, S, FA consider the model formulations where nucleotides catalyse carbon fixation and energy, sugar or fatty acid production respectively.


#

Running the Fig3 script results in Figure 3. Fig2 sweeps through different Km^N parameter values in the different 'intFunctionNucCat' functions to integrate the equations present in the  'partFunctionNucCat' functions. There is one version of these functions for nucleotide catalysis the synthesis pathways of each species (FA,AA,S and E), as indicated by the name of the species in the end of the function name.

This is done using K_C2_N = 10^0.9 and cost = 0 for Figure 3a, K_C2_N = 10^3 and cost = 0 for Figure 3b, and K_C2_N = 10^3 and cost = 1 for Figure 3c.


## Figure 4


Running the Fig4 script results in Figure 4. This is done much in the same way as Figure 3, but the integration functions have concentration of nucleotides as an output instead.



## Figure 5 
Running the Fig5 script results in Figure 5. Fig5 sweeps through a range of parameters for the rate of nucleotide catalysis of carbon fixation and the rate of nucleotide autocatalysis, using the 'intFunctionNucCatC2_parsweep' function to integrate the equations in 'partFunctionNucCatC2_parsweep'. It then uses the results of the simulations to plot two heatmaps. 
