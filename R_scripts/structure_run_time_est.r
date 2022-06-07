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
