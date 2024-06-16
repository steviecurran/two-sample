# two-sample
Two sample test 

Python code to give the confidence interval for two sample test.

Uses z and t-distributions, without the need of look-up tables

Data can be extracted any column in an ascii for csv files.

E.g. Data wtih the two samples to compare? 3.15_apples.dat

For NY -  n = 10, mean = 3.941, SD = 0.184 [sample]

For LA -  n = 8, mean = 3.245, SD = 0.268 [sample]

Level of confidence [e.g. 95, 99, 99.9% - z = 3 sigma is 99.75]? 95

At least one sample size < 30 so using t-value

Variance ratio is 0.47,

Assume equal sample variances (for ratio between 0.5 and 2 can assume equal) [y/n]? y

For 95.00% confidence (16 DoFs), t-value is 2.120, giving mean diff of 0.70 +/- 0.23 (0.47 to 0.92)

-----------------------------------------------------------------------------------------------------

For one sample version see https://github.com/steviecurran/one-sample

For 99.90% confidence (16 DoFs), t-value is 4.006, giving mean diff of 0.70 +/- 0.43 (0.27 to 1.12)




