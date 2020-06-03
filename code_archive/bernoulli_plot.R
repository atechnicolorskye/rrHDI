## Reproduce results for Bernoulli 2013: Ridge Projection for hdi
rm(list=ls())

load("out/sim_hdi_sims10_n100.rda")

M1 = Results[Results$Sig == "id", ] 
M2 = Results[Results$Sig == "equi.corr", ] 

pdf(file='fig/bernoulli_2.pdf')
par(mfrow=c(2,2))
p_single_M1 = plot((M1$fwer_single+1e-8), M1$pwr_single, log="x", xlim=c(1e-8, 0.1),
                   main = "single testing (M1)", ylab = "average power", xlab = "average type I error")
abline(v=0.05,lty=2,col="black")
p_single_M2 = plot((M2$fwer_single+1e-8), M2$pwr_single, log="x", xlim=c(1e-8, 0.1),
                   main = "single testing (M2)", ylab = "average power", xlab = "average type I error")
abline(v=0.05,lty=2,col="black")
p_multi_M1 = plot((M1$fwer_multi+1e-8), M1$pwr_multi, log="x", xlim=c(1e-8, 0.1),
                  main = "multiple testing (M1)", ylab = "average power", xlab = "FWER")
abline(v=0.05,lty=2,col="black")
p_multi_M2 = plot((M2$fwer_multi+1e-8), M2$pwr_multi, log="x", xlim=c(1e-8, 0.7),
                  main = "multiple testing (M2)", ylab = "average power", xlab = "FWER")
abline(v=0.05,lty=2,col="black")
dev.off()
