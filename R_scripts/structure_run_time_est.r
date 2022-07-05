#estimated run time from Raj (faststructure paper) fig 4 in sqrt sqrtregression

sqrt_prob_size = c(500, 1000, 1500, 2000)
sqrt_run_time = c(5,10,15,20)
#slope is 0.01

m = 0.01
b = 0

k = c(1:15)
n = 102 # number samples
l = 40857 # number loci

#gives estimated run time in days
(( sqrt(k*n*l) * m + b)^2 ) / 60 / 24


[1] 0.2894037 0.5788075 0.8682112 1.1576150 1.4470188 1.7364225 2.0258263
 [8] 2.3152300 2.6046338 2.8940375 3.1834413 3.4728450 3.7622487 4.0516525
[15] 4.3410563
