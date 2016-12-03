# Point-Estimate-Methods

Point Estimate Methods can be used to approximate the nth order statistical moments of a random variable Y, where Y = f(X) and X is the input random variable to a function. They do so by calculating points on X with a certain "significance", and attaching weights to those points. These points and weights are then used in the function to output points, y_i = f(x_i), and the n-th order statistical moments are calculated by using these y_i and the weights.

Point Estimate Methods are significant in places where Monte-Carlo simulations require a long period of time, because instead of using a large number of points to perform a simulation, a PEM method can be used to speed things up since it only requires 2m or 2m+1 points, where m is the number of dimensions of the random variable.

The File "main_for_pem" illustrates this idea by calculating the first three central moments of 3 different distributions, and comparing them to analytic/Monte-Carlo solutions, by presenting a table. It analysizes the following PEM schemes:

1. Rosenblueth's PEM
2. Hong's 2-point and 3-point estimate schemes.
3. Zhao's 7 point estimate scheme.

Results are also presented graphically to visualize the points and their weights chosen from the distribution. Zhao's PEM utilizes the inverse Rosenblatt transform, which is calculated in the PEM function.

#Sampling

Two different sampling methods have been implented in the "samplers" function file. The metropolis-hastings sampler, based on random walks is implemented, and a lating hyper cube sampler has also been implemented. It is easy to notice that LHS provides a more balanced sample (with repsect to the pdf), and metropolis-hastings requires burn-in iterations and doesn't necessarily provide a very balanced sample size, especially for a lower number of samples.

These two functions (PEM and Sampling) were implemented as part of my work at the KIT institute of Applied Informaticks, for analyses of Stochastic Optimal Power Flow methods. The team webpage and more information about them can be found here:
https://www.iai.kit.edu/english/1399.php
