This project contains a very simple simulation to explore two implementations
of Beyesian mixed modeling. The simulation creates data with a linear
relationship between one covariate and a response variable with an intercept
that varies by site nested within species.

The Stan code is a literal implementation of the simulation code with just a
few small changes:
<ul>
<li>It uses a <tt>cauchy(0,5)</tt> prior on the residual variance.</li>
<li>It uses a <tt>gamma(1,1)</tt> prior on the random effect variances.</li>
The <tt>stan_lmer()</tt> model is the the direct analog of the simulation and
the Stan code.
</ul>

<tt>results.txt</tt> contains the results of simulations under several
different prior specifications for the random effects. You'll notice
that the bias, root mean squared error, and coverage (symmetric 80%
credible intervals) are very similar for the intercept
(<tt>beta_0</tt>), the regression coefficient (<tt>beta_1</tt>), nad
the residual standard deviatiion (<tt>sigma</tt>). The random effect
standard deviations (<tt>sigma_species</tt> and
<tt>sigma_species_site</tt>) are fairly different from one another and
the Stan code has a smaller bias and root mean squared error and it
has better coverage properties.