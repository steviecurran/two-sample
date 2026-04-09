## A/B Testing Toolkit: Comparing Two Means

Many teams run experiments but struggle with statistical validity and decision-making under uncertainty. This toolkit provides the means to 
determine whether an observed difference between groups is statistically significant and practically meaningful, as commonly required in product experiments and business decision-making.

The notebook, *CI-2_means.ipynb*, compares two groups using confidence intervals and statistical testing, without lengthy manual calculations and distribution tables which can lead to errors.

The notebook is interactive and includes options to:

- Input the full dataset  in .csv and .dat format or input the means and standard deviations directly

- For the full data option, plot a histogram showing the distributions of the two classes 

- Adjust the confidence levels from the default 95%, offering more stringent testing

- Run a one or two tailed test

- Change the t- to z-statistic threshold from the default n = 30

- Switch between equal and inequal variances
 
**Example 1: Raw data of a small sample**

The file *Mg_levels.dat* contains the levels of magnesium in a sample of people before and after taking a supplement. We wish to test the null hypothesis that the supplement does not increase the magnesium levels of the patients

![](https://raw.githubusercontent.com/steviecurran/two-sample/refs/heads/main/Mg_1.png)

![](https://raw.githubusercontent.com/steviecurran/two-sample/refs/heads/main/Mg_histo.png
)

![](https://raw.githubusercontent.com/steviecurran/two-sample/refs/heads/main/Mg_results.png)

So at 95% confidence we can reject the null hypothesis that the supplements do not increase the patient's magnesium levels.

Using the condifence level dropdown it is straightforward to show that the null hypothesis cannot be rejected with 99% confidence.

**Example 2: Large sample test**

We can use the medical data below to demonstrate the toolkit on a large sample where we input the summary data only.

![](https://raw.githubusercontent.com/steviecurran/two-sample/refs/heads/main/medical_large.png)

For example, for the *systolic blood pressure*

![](https://raw.githubusercontent.com/steviecurran/two-sample/refs/heads/main/bp_results.png)

Here *Men* have been entered as Sample 1 and *Women* as Sample 2.

So although the means are only about 1/3 of standard deviation apart, we are 95% confident that men have  higher systolic blood pressure.  This holds above 99.999% confidence.
 
