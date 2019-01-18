<h1>ssev</h1>

<p>The goal of ssev is to compute optimal sample sizes for two-stage RCT designs when incorporating the size of the population. The basic idea is simple: given an estimate of the effect size of the difference between two groups, and an estimate of the population size N, we can compute the optimal size of an RCT that maximizes the expected outcome over both the RCT and the resulting guideline that is used for the remainder of the population. The package computes the associated sample size for a given significance level when comparing two means (with equal or unequal variances assumed) or when comparing two proportions. </p>

<h2>Example</h2>

<p>The following code computes the sample size required for a population of \(N=10^5\) units when assuming equal variances and a cohen&#39;s d of .8:</p>

<pre><code class="r example">compute_sample_size(means=c(0,.8), sds=1, N=10^5)
</code></pre>
