# Meeting Notes:

## Meeting 25

### Meeting Date: 18/10/18

### Who was present:
* David Dobbie
* Paul Teal

### Meeting Comments/Notes:
* finishing final report
* will practice presentation next week w/ Paul

### Action points for next meeting:
* full report done
* presentation draft done

### Action points achieved since last meeting:
* wrote most of report
* implemented a compromise of ILT+ so that it may be indicative

### Action points yet to be achieved since last meeting:
* almost finished report



## Meeting 24

### Meeting Date: 11/10/18

### Who was present:
* David Dobbie
* Paul Teal

### Meeting Comments/Notes:
* problem with ILT+ comes from the moment estimator

### Action points for next meeting:
* finish report by then
* validate ILT+ as much as possible

### Action points achieved since last meeting:
* checked ILT+ validation once again

### Action points yet to be achieved since last meeting:
* written conclusion of report 1st draft




## Meeting 23

### Meeting Date: 4/10/18

### Who was present:
* David Dobbie
* Paul Teal

### Meeting Comments/Notes:
* got comments back for report sans conclusion

### Action points for next meeting:
* write conclusion (3 pages)
* check ILT+ validation, at least near more near to a degree

### Action points achieved since last meeting:
* written most of report, gone through first draft process

### Action points yet to be achieved since last meeting:
* written conclusion of report 1st draft


## Meeting 22

### Meeting Date: 27/09/18

### Who was present:
* David Dobbie
* Paul Teal

### Meeting Comments/Notes:
* setting up test environment, polishing test environment

### Action points for next meeting:
* write draft 2nd half

### Action points achieved since last meeting:
* write draft 1st half
* sorted out comparison set up for validating recreated models

### Action points yet to be achieved since last meeting:
* none





## Meeting 21

### Meeting Date: 20/09/18

### Who was present:
* David Dobbie
* Paul Teal

### Meeting Comments/Notes:
* validating ILT+ and ILT. Can get to within R2 = 0.9 for ILT+. Can get close but not completely on.
* can add case by case comparisons with known ILT and ILT+ in Gruber 2013. Only require BFF to be accurate


### Action points for next meeting:
* start drafting 2nd half of report
* set up comparison cases for models presented in Gruber 2013. We can make a comparison there.


### Action points achieved since last meeting:
* written most of first half first draft, aim to hand in on Friday 21/09
* set up validation environment for ILT and ILT+ methodologies

### Action points yet to be achieved since last meeting:
* ILT+ validation yet to be fully 100% completed


## Meeting 20

### Meeting Date: 06/09/18

### Who was present:
* David Dobbie
* Paul Teal

### Meeting Comments/Notes:
* possible chance of taking practical T2 relaxation data to further improve on 489 report
* started validation of the ILT and ILT+ techniques for the test environment
* started integration with preliminary report, moving onto writing up design chapter

### Action points for next meeting:
* double-triple check ILT+ technique, make evaluation as fair as possible
* Start writing report, 
  * hand in  Design sections draft 1
* send email to investigate possible experimental data from a 


### Action points achieved since last meeting:
* in progress to complete checking ILT+ and writing design chapter

### Action points yet to be achieved since last meeting:
* most of drafting for the first half of the report has been finished





## Meeting 19

### Meeting Date: 30/08/18

### Who was present:
* David Dobbie
* Paul Teal

### Meeting Comments/Notes:
* look into creating Parzan window (mixed Gaussian) rather than simple Gaussian for each data point
* start writing the report, combined portions of prelim into full
* keep in mind feedback for the prelim report when writing

### Action points for next meeting:
* double-triple check ILT+ technique, make evaluation as fair as possible
* Start writing report, 
  * hand in Intro, Background, Design sections draft 1 (mid-tri wk 2)
  * hand Implementation, Evaluation, Conclusion sections draft 1 (wk7)


### Action points achieved since last meeting:
* set up the BBF estimator for:
  4. ILT+ 2013 technique
  5. Tapered area estimator 2013



### Action points yet to be achieved since last meeting:
* none








## Meeting 18

### Meeting Date: 23/08/18

### Who was present:
* David Dobbie
* Paul Teal

### Meeting Comments/Notes:
* gone over testing system for the BFF estimator
* Could create a Parzan window  for combining the prior - a Gaussian mixture

### Action points for next meeting:
* set up the BBF estimator for:
  4. ILT+ 2013 technique
  5. Tapered area estimator 2013

### Action points achieved since last meeting:

* set up the BFF estimator for:
  1. Bayesian with steep cut off for fluid volume
  2. Bayesian with tapered cut off for fluid volume
  3. Classic 2002 Venk technique
* set up testing environment


### Action points yet to be achieved since last meeting:
* set up the BBF estimator for:
  4. ILT+ 2013 technique
  5. Tapered area estimator 2013


## Meeting 17

### Meeting Date: 16/08/18

### Who was present:
* David Dobbie
* Paul Teal

### Meeting Comments/Notes:
* tested how Gaussian the prior data points are.
* set up a One-sample Kolmogorov-Smirnov test for each porosity bin in f(T2) - finding that a collection of points have a Gaussian distribution
* looked into testing the log-normal closeness of the density function prior data points. This process is still pending

### Action points for next meeting:
* set up the BBF estimator for:
  1. Bayesian with steep cut off for fluid volume
  2. Bayesian with tapered cut off for fluid volume
  3. Classic 2002 Venk technique
  4. ILT+ 2013 technique
  5. Tapered area estimator 2013
* set up testing environment for the system

### Action points achieved since last meeting:
* tested how Gaussian each data point of the prior density function is
* finished off testing different estimate of covariance matrices

### Action points yet to be achieved since last meeting:






## Meeting 16

### Meeting Date: 09/08/18

### Who was present:
* David Dobbie
* Paul Teal

### Meeting Comments/Notes:
* explored the performance of different estimators for the covariance matrix
* we see that there are diminishing returns at certain estimation of alpha (how much we trust the prior)
* the non diagonal covariance estimator indicates significant performance improvements and less sensitivity for the choice of the alpha (it is around 1, which fits with the hypothesis that there should be equal weighting for that)
* for an analysis of the cut-off time for bound fluid: it would be best served to choose a more arbitrary point so that there can be more useful discussion of the model's effectiveness.


### Action points for next meeting (or two):
* test how Gaussian each data point of the prior density function is
  * how true is the model's assumption that the density function is is normal for each point
* create an experiment of the zero uniform offset on the estimated covariance matrix
* set up the BBF estimator for:
  1. Bayesian with steep cut off for fluid volume
  2. Bayesian with tapered cut off for fluid volume
  3. Classic 2002 Venk technique
  4. ILT+ 2013 technique
  5. Tapered area estimator 2013

### Action points achieved since last meeting:
* explored the use of different covariance matrices for the Bayesian estimator.
  * recorded results for 1.3, 1.4, 1.5 version of the Bayesian estimator
* developed a simple alpha predicted with the 2 norm of the covariance matrix


### Action points yet to be achieved since last meeting:
* have not made a system that trained an alpha for the estimator - far too high complexity for something that can be simply done with a more suitable covariance matrix






## Meeting 15

### Meeting Date: 2/08/18

### Who was present:
* Paul Teal
* David Dobbie

### Meeting Comments/Notes:
* have presented the preliminary presentation. Good feedback given for creating an accessible presentation
* we have a framework where experimental data can be used to form the estimator
  * this uses leave one out cross validation to give an empirical best alpha
* the path for the project is to use different priors formed by experimental data such as
  * covariance of each data point
  * assuming smoothness, the cross-correlation between each data point
* also given feedback on the preliminary report
  * give discussion of results and experiments
  * discuss how algorithms had their recreation validated
  * first letter capital for figures and sections cross referencing (Figure 1, Section 2.1)
  


### Action points for next meeting:
* create a method for finding the best alpha value. Try two methods:
  * create a validation on the training process and choose the optimal empirical alpha value from that
  * try to estimate alpha as we already know C_f and C_n - see if this matches

### Action points achieved since last meeting:
* presented the preliminary presentation
* implemented a cross validation system using experimental data to evaluate performance for different alpha values - do not need to normalise data before using it
  * this implementation uses a uniform diagonal covariance for the density function

### Action points yet to be achieved since last meeting:
* test how normal the density function experimental data is
* create weighted covariance C_f dependent on previous prior data



## Meeting 14

### Meeting Date: 26/07/18
### Who was present:
* Paul Teal
* David Dobbie

### Meeting Comments/Notes:
* was sick for a significant part of this week, could not present
* dealing with a normalisation issue of the experimental data. Normalisation directly assumes the porosity itself.
  * will construct prior with normalised data first so bound fluid fraction can be computed
* could perhaps utilise hyper priors for density estimation

### Action points for next meeting:
* normalise and construct the prior
* preliminary presentation in wk 3
* test how normal original data is
* estimate covariance, use it to make prior

### Action points achieved since last meeting:
* loaded experimental data onto a script to construct the prior

### Action points yet to be achieved since last meeting:
* full construction of a Gaussian prior
* preliminary presentation


## Meeting 13

### Meeting Date: 18/07/18

### Who was present:
* Paul Teal
* David Dobbie

### Meeting Comments/Notes:
* constructed the Bayesian technique with a zero mean assume prior density function
* could perhaps rework the density function distribution to be log normal to work around the non-negative constraint on that
* the results gathered use a realistic prior with using the true density function as the mean. This is the best possible circumstance
* will present in front of an ENGR 489 lecture slot for the mid way presentation
* have not received preliminary report results yet

### Action points for next meeting:
* construct a prior with high quality experimental data, use for Bayesian technique (assume independence of each data point). This will be compared with synthesised measurement data created from it
* test how Gaussian (or not Gaussian) the experimental density functions are (normality test)
* look into improving the prior of experimental data by creating a realistic covariance (for correlation between points) (use autocorrelation method and covariance method)

### Action points achieved since last meeting:
* implemented and evaluated Bayesian technique with a zero mean prior density function

### Action points yet to be achieved since last meeting:
* create trimester 2 presentation






## Meeting 12

### Meeting Date:  28/06/18

### Who was present:
* Paul Teal
* David Dobbie

### Meeting Comments/Notes:
* Paul will not be able to meet for next two Thursdays (5/7/18 and 12/7/18)
* Was unable to meet up due to exams - compensated expectation as a result
* discussion about MSE in evaluating estimators

### Action points for next meeting:
* Start Paul's Bayesian technique - construct and evaluate it
* Sort out T2 week one presentation

### Action points achieved since last meeting:
* Mostly finished Venkataramnan et al Porosity Estimation (April 2015)
  * understanding of how the method works as well as its limitations

### Action points yet to be achieved since last meeting:
* none


## Meeting 11

### Meeting Date: 7/06/18

### Who was present:
* Paul Teal
* David Dobbie

### Meeting Comments/Notes:
* generate the bias error from individual deltas = f(T2) for all of T2

### Action points for next meeting:
* Finished Venkataramnan et al Porosity Estimation (April 2015)

### Action points achieved since last meeting:
* progress on Venkataramnan et al Porosity Estimation (April 2015)
* final tweaking of preliminary report completed

### Action points yet to be achieved since last meeting:
* none



## Meeting 10

### Meeting Date: 31/05/18

### Who was present:
* Paul Teal
* David Dobbie

### Meeting Comments/Notes:

### Action points for next meeting:
* Tweak preliminary report, more detail on the future plan and proposed technique
  * rearrange for work done
  * add more questions to the marker
* Start Venkataramnan et al Porosity Estimation (April 2015)

### Action points achieved since last meeting:
* completed full draft 

### Action points yet to be achieved since last meeting:
* none


## Meeting 9

### Meeting Date: 24/05/18

### Who was present:
* Paul Teal
* Robin Dykstra
* David Dobbie

### Meeting Comments/Notes:
* first half draft returned, needs restructuring and rewriting
* requires specific referencing


### Action points for next meeting:
* finished first full draft on Monday - (fully use feedback given to improve it)

### Action points achieved since last meeting:
* Finished T2 ILT+ Gruber 2013
* Fixed figures in report
* finished first draft of other half of report, ready to combine together into effective script

### Action points yet to be achieved since last meeting:
* None




## Meeting 8

### Meeting Date: 17/05/18

### Who was present:
* Paul Teal
* David Dobbie

### Meeting Comments/Notes:
* acquired Hurlimann chapter on oil well logging
* format report so that:
   * smallest text is the size of a sub caption
   * remove titles on figures
   * refer to subplots by letters, allows for elegant discussion
   * implement hierarachal writing structure - each chapter and section have a description at the start describing what they will entail.
* Robin will come into next meeting to help out with writing the background and motivation on the NMR material.

### Action points for next meeting:
* write up introduction, future work and ILT+ method of preliminary report
* finish debugging Gruber 2013, paper 4
* change formatting of report, make it better

### Action points achieved since last meeting:
* Wrote 80% of background of the preliminary report
* Created script evaluating tapered area and moment estimations with the true values
* created script framework implementing ILT+ and ILT together


### Action points yet to be achieved since last meeting:
* ILT+ not working yet, need to fix the weighting matrix - make it work as intended


## Meeting 7

### Meeting Date: 10/05/18

### Who was present:
* Paul Teal
* David Dobbie

### Meeting Comments/Notes:
* emphasis on images and diagrams for communicating concepts
* place evaluation in terms of mean squared error - easy to keep track of
* ask for Haliburtion Oil Well Logging book from Robin for useful images
* bring USB for Martin Hurlimann chapter from Paul with useful images



### Action points for next meeting:
* Half of preliminary report complete - the introduction and background survey
* Finish Gruber Paper 4
* Gather accuracy of tapered areas estimations (for report)


### Action points achieved since last meeting:
* Finished Gruber Paper 3 (Tapered Areas)


### Action points Yet to be achieved since last meeting:
* none



## Meeting 6

### Meeting Date: 03/05/2018

### Who was present:
* Paul Teal
* David Dobbie

### Meeting Comments/Notes:
* Given CPMG measurement data - not too sure what method was used to derive the T2 distributions from it
* Paper 3 gives tapered areas and moments to establish the constraints of paper 4 - Gruber  T2 distribution
* A major part of this project is to establish Tc - the fraction cut off point between the bound and free fluids in the sample. Paper 3 discusses this
   * good to have kernel made specially to differentiate Tc cut off point
* can ignore the incomplete polarisation (G(t)) modelling in S2.2 of paper 3, want to constrain scope to fully polarised measurements (M(t)).
* do not need to have complete readings of background material due to constraints on project execution. Very direct implementation

### Action points for next meeting:
* Complete Gruber 2013 paper 3
* Start Gruber 2013 paper 4 (optimistic goal to complete it)

### Action points achieved since last meeting:
* written up summaries of venk 2002 and venk 2010
* started gruber 2013 on linear functionals.

### Action points Yet to be achieved since last meeting:
* none



## Meeting 5

### Meeting Date: 
26/04/18

### Who was present:
* David Dobbie
* Paul Teal

### Meeting Comments/Notes:
* Stopping work on Mellin Transform - got satisfactory results for omega in [0,1]
* making sure documentation for past papers are good - able to revise upon on final report
* moving onto twin papers by Gruber 2013 that builds on this
 
### Action points for next meeting:
* Write up summary of venk 2002 and venk 2010
* Start Gruber 2013 - Estimation of Petrophysical and fluid properties using integral transforms in nuclear magnetic resonance

### Action points achieved since last meeting:
* Finished implementing mellin transform venk 2010, written documentation on the process.

### Action points Yet to be achieved since last meeting:




## Meeting 4

### Meeting Date: 
19/04/18

### Who was present:
* Paul Teal
* David Dobbie 

### Meeting Comments/Notes:
* look into using logarithms to get around computational errors and problems around float computation problems

### Action points for next meeting:
* Bug check implementation of Mellin transform Venk.2010, have pictures and results of bug testing (bug plots of G(\omega) and \delta_i) (before meeting)


### Action points achieved since last meeting:
* progress on implementation of Venk. 2010



### Action points yet to be achieved since last meeting:
* Not yet completed implementation of Venk. 2010 Mellin transform








## Meeting 3

### Meeting Date: 
13/04/18

### Who was present:
* Paul Teal
* David Dobbie 

### Meeting Comments/Notes:
* Would in future look into implementing the flint code, play around with it to make it work appropriately
* important to recreate the paper's properly to make them comparable
* Can utilise "engauge" to convert image data to data points to recreate the paper's made
* Further on with Bayes, about constructing the Priori for the Bayesian analysis

### Action points for next meeting:
* Work onto implementing Mellin Transform of CPMG data (Journal of Magnetic Resonance 206 (2010) 20-31) (note that this is more of a stepping stone paper)


### Action points achieved since last meeting:
* finished 1D implementation of Venk. 2002
* Progress on 2D implementation of Venk. 2002


### Action points yet to be achieved since last meeting:
* none




## Meeting 2

### Meeting Date: 
06/04/18
### Who was present:
* Paul Teal
* David Dobbie

### Meeting Comments/Notes:
* got code to guide implementation of 1D example - aiming to implement K1 kernel properly, etc.
* Figured out aim (ambitious) to code up competing algorithms and understand by preliminary report hand in.
* given series of papers for next few weeks of work, get a comprehensive picture of the project.
* Discussion on introductory Bayes, will be focus later on in the project.

### Action points for next meeting:
* Finish 1D example of extracting exponential time constants in paper 1 venkataramanan-2002-tensor-fredholm
* Get strong understand of 2D example of it, implement it

### Action points achieved since last meeting:
* Read paper (Vol. 50, no 5. May 2002 IEEE transactions on Signal Processing (1053-587x(02)03282-8) venkataramanan-2002-tensor-fredholm.pdf
* Worked on implementing the algorithm in it, not complete

### Action points Yet to be achieved since last meeting:
* Implement 1D example fully in paper venkataramanan-2002-tensor-fredholm




## Meeting 1

### Meeting Date:
 20/03/18
### Who was present:
* Paul Teal
* Robin Dykstra
* David Dobbie

### Meeting Comments/Notes:
-	Going over project knowledge base
-	With the project, techniques required to compensate for the very poor SNR in underground oil detection with how NMR used
- A magnetic field when constant results in a constant spin characteristic of the protons itself (\omega = \beta\gamma)
- When sending a pulse to redirect the spin (90 degrees for example), we have the regrow to normal spin T1 and the much faster decay in T2.
- We have T2* for a non-constant magnetic field, dealing with decay caused by a non-uniform magnetic field since the 2D vectors drift until they oppose.
- To get around this faster decay, we send a 180 degree pulse of RF radiation to reorient it so that it acts in the envelope of the T2 decay curve – we get discretised measurements of T2. (these are called echos)
-	When a molecule is closer to the surface, for a fluid (bound) rather than deeper inside (free), it will decay faster
-	We can get different t2s for different protons and their different physical positioning
-	If we take several measurements of T2 exponential decays we can make a statistical reference of this data and see the prominence of both (we have sum of exponentials with different modifiers on T2)
-	Looking at either T2 or T1 combined gives a 2.5D representation of it, this is in the paper Vol. 50, no 5. May 2002 IEEE transactions on Signal Processing (1053-587x(02)03282-8)
- This paper utilised solving for the distribution from discretised data with the Kroncher product
- 	Note that SNR is defined differently in this paper than its conventional meaning.
-	Initially, will try to get the distributions sorted for a very simplified version of this problem, look up the Bulter-Reed-Dawson algorithm
- WILL NEED to consider a priori that can allow the data to converge so it can be used properly.
-	Given book: Principles of Magnetic Resonance Imaging (A Signal Processing Perspective), look at eq 3.69 (the Bloch equation), it is the effect of the applied magnetic field used for MRI
-	Project will be done in MATLAB

### Action Points For Next Meeting:
- Read the (1053-587x(02)03282-8) IEEE paper, get understanding of what is happening
- Implement Butler-Reed-Dawson algorithm and a simplified version of the method in the paper  with 1 dimension for T2
