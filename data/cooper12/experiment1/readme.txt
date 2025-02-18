This folder contains data from Experiment 1 of Cooper, Wutke & Davellar (2012)
"Differential contributions of set-shifting and monitoring to dual-task
interference", Quarterly Journal of Experimental Psychology, 65 (3), 587-612.

The data are arranged as follows:

auxiliary_data.tab is a tab-delimited data file that contains the following
columns:
 * partnr: The unique participant identifier
 * group: Always 1 for this experiment
 * sex: coded as 1 for male and 2 for female
 * dob: Year of birth - used to calculate age in 2008 when the study was done
 * educ: Level of education (from 1 = high school to 5 = PhD)
 * DS_n: Number of trials on the digit switching task
 * DS_sw: Number of switch errors
 * DS_miss: Number of trials with no response
 * 2B_n: Number of trials on the 2-back task
 * 2B_hit: Number of hits on the 2-back task
 * 2B_cr: Number of correct rejections on the 2-back task
 * 2B_fp: Number of false positives on the 2-hack task
 * 2B_miss: Number of misses on the 2-back task
 * GN_n: Number of trials on the go/no-go task
 * GN_hit: Number of hits on the go/no-go task
 * GN_cr: Number of correct rejections on the go/no-go task
 * GN_fp: Number of false positives on the go/no-go task
 * GN_miss: Number of misses on the go/no-go task
Data from each participant is on a separate line.

Data from the RNG task is contained in set of files with names of the form
rng_XXXX_Y.txt where XXXX is the participant number (ranging from 0000 to
0036) and Y encodes the condition (0 = control; 1 = with digit switching; 2 =
with 2-back; 3 = with go/no-go). Each of these files consists of 100 lines
with successive lines encoding details of successive responses. Each line
should contain two numbers: an integer representing the respone (0 = 'A' up to
9 = 'J') and a real number that gives the time between the previous response
and this response (in seconds). A response of '-1' indicates no response
within the timeout perioud (5 seconds).

Richard P. Cooper
Mon Sep 14 14:15:29 2015


