## Reproduce results for Bernoulli 2013: Ridge Projection for hdi
rm(list=ls())

load("out/sim_blpr_sims1000_b25.rda")

cp = Results$covprob
intlen = Results$avgintlen

pdf(file='fig/blpr.pdf', width = 8, height = 5)
par(mfrow=c(1,2))
p_cp = plot(1:25, cp, ylim=c(0.8, 1.0), main = "Coverage probability", ylab = "coverage probability",
            xlab = "index")
p_il = plot(1:25, intlen, ylim=c(0.0, 0.5), main = "Interval length", ylab = "interval length",
            xlab = "index")
dev.off()