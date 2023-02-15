rm(list = ls())

library('svglite')
source('../functions/ggsave2x.R')
source('../functions/dfcombos.R')


res <- 0.01

T_opt0 <- 18
T_min0 <- 0
T_max0 <- 25
qhat_opt0 <- 1.5
qhat0 <- function(temp_C0) {
  qhat0 <- qhat_opt0 * ((temp_C0 - T_max0) * (temp_C0 - T_min0)^2) / 
    ((T_opt0 - T_min0) * ((T_opt0 - T_min0) * (temp_C0 - T_opt0) - 
                            (T_opt0 - T_max0) * (T_opt0 + T_min0 - 2*temp_C0)))
  return(qhat0)
}
temp_C0 <- seq(0,25,res)


res <- 0.01

T_opt <- 18
T_min <-10
T_max <- 25
qhat_opt <- 3.6
qhat <- function(temp_C) {
  qhat <- qhat_opt * ((temp_C - T_max) * (temp_C - T_min)^2) / 
    ((T_opt - T_min) * ((T_opt - T_min) * (temp_C - T_opt) - 
                          (T_opt - T_max) * (T_opt + T_min - 2*temp_C)))
  return(qhat)
}
temp_C <- seq(10,25,res)

res <- 0.01

T_opt_old <- 18
T_min_old <-0
T_max_old <- 25
qhat_opt_old <- 3.6
qhat_old <- function(temp_C_old) {
  qhat_old <- qhat_opt_old * ((temp_C_old - T_max_old) * (temp_C_old - T_min_old)^2) / 
    ((T_opt_old - T_min_old) * ((T_opt_old - T_min_old) * (temp_C_old - T_opt_old) - 
                          (T_opt_old - T_max_old) * (T_opt_old + T_min_old - 2*temp_C_old)))
  return(qhat_old)
}
temp_C_old <- seq(0,25,res)


T_opt2 <- 28
T_min2 <- 10
T_max2 <- 38
qhat_opt2 <- 5.6
qhat2 <- function(temp_C2) {
  qhat2 <- qhat_opt2 * ((temp_C2 - T_max2) * (temp_C2 - T_min2)^2) / 
    ((T_opt2 - T_min2) * ((T_opt2 - T_min2) * (temp_C2 - T_opt2) - 
                            (T_opt2 - T_max2) * (T_opt2 + T_min2 - 2*temp_C2)))
  return(qhat2)
}
temp_C2 <- seq(10,38,res)

T_opt3 <- 36
T_min3<- 15
T_max3 <- 45
qhat_opt3 <- 7.2
qhat3 <- function(temp_C3) {
  qhat3 <- qhat_opt3 * ((temp_C3 - T_max3) * (temp_C3 - T_min3)^2) / 
    ((T_opt3 - T_min3) * ((T_opt3 - T_min3) * (temp_C3 - T_opt3) - 
                            (T_opt3 - T_max3) * (T_opt3 + T_min3 - 2*temp_C3)))
  return(qhat3)
}
temp_C3 <- seq(15,45,res)

T_opt4 <- 43.75
T_min4<- 26.25
T_max4 <- 51.25
qhat_opt4 <- 8
qhat4 <- function(temp_C4) {
  qhat4 <- qhat_opt4 * ((temp_C4 - T_max4) * (temp_C4 - T_min4)^2) / 
    ((T_opt4 - T_min4) * ((T_opt4 - T_min4) * (temp_C4 - T_opt4) - 
                            (T_opt4 - T_max4) * (T_opt4 + T_min4 - 2*temp_C4)))
  return(qhat4)
}
temp_C4 <- seq(26.25,51.25,res)

T_opt5 <- 55
T_min5<- 30
T_max5 <- 60
qhat_opt5 <- 8
qhat5 <- function(temp_C5) {
  qhat5 <- qhat_opt5 * ((temp_C5 - T_max5) * (temp_C5 - T_min5)^2) / 
    ((T_opt5 - T_min5) * ((T_opt5 - T_min5) * (temp_C5 - T_opt5) - 
                            (T_opt5 - T_max5) * (T_opt5 + T_min5 - 2*temp_C5)))
  return(qhat5)
}
temp_C5 <- seq(30,60,res)



svglite("../figures/appendix_qhat_dynamics.svg", width = 16/2.54, height = 8/2.54)

par(mfrow=c(1,2))
par(mar=c(4,5,1,1))
par(cex.axis=1, cex.lab=1, cex.main=0.6, cex.sub=0.6)

plot(temp_C, qhat(temp_C),xlim = c(0,60), ylim=c(0, 10), type ="l", xlab = '',
     ylab = expression(q['max,opt']~(g[COD-S]~g[COD-B]^'-1'~d^'-1')), col = "blue", lwd = 1, lty = 1)
lines(temp_C0, qhat0(temp_C0), col = "black", lwd = 1, lty = 2)
lines(temp_C2, qhat2(temp_C2), col = "orange", lwd = 1, lty = 2)
text(18, 4.0, "18, 3.6", cex = 0.6)
text(18, 1.7, "18, 1.5", cex = 0.6)
text(28, 6, "28, 5.6", cex = 0.6)

legend(0,10, legend = c('m0', 'm1_mod', 'm2'), col = c('black', 'blue', 'orange'), lty = c(2, 1, 2), box.lwd = '', cex = 0.6)
title('this work')

plot(temp_C_old, qhat_old(temp_C_old),xlim = c(0,60), ylim=c(0, 10), type ="l", xlab = '',
     ylab = expression(q['max,opt']~(g[COD-S]~g[COD-B]^'-1'~d^'-1')), col = "blue", lwd = 1,lty = 1, cex.axis = 1, cex.lab = 1)
lines(temp_C2, qhat2(temp_C2), col = "orange", lwd = 1, lty = 2)
lines(temp_C3, qhat3(temp_C3), col = "palegreen4", lwd = 1,lty = 1)
lines(temp_C4, qhat4(temp_C4), col = "purple", lwd = 1, lty = 2)
lines(temp_C5, qhat5(temp_C5), col = "red", lwd = 1, lty = 1)

text(18, 4.0, "18, 3.6", cex = 0.6)
text(28, 6.0, "28, 5.6", cex = 0.6)
text(36, 7.6, "36, 7.2", cex = 0.6)
text(43.75, 8.4, "44, 8", cex = 0.6)
text(55, 8.4, "55, 8", cex = 0.6)

legend(0,10, legend = c('m1','m2','m3','m4','m5'), col = c('blue', 'orange','palegreen4','purple','red'), lty = c(1, 2, 1,2,1), box.lwd = '', cex = 0.6)
title('original')

dev.off()

png('../figures/appendix_qhat_dynamics.png',  width = 16/2.54, height = 8/2.54, units = 'in', res = 600)

par(mfrow=c(1,2))
par(mar=c(4,5,1,1))
par(cex.axis=1, cex.lab=1, cex.main=0.6, cex.sub=0.6)

plot(temp_C, qhat(temp_C),xlim = c(0,60), ylim=c(0, 10), type ="l", xlab = '',
     ylab = expression(q['max,opt']~(g[COD-S]~g[COD-B]^'-1'~d^'-1')), col = "blue", lwd = 1, lty = 1)
lines(temp_C0, qhat0(temp_C0), col = "black", lwd = 1, lty = 2)
lines(temp_C2, qhat2(temp_C2), col = "orange", lwd = 1, lty = 2)
text(18, 4.0, "18, 3.6", cex = 0.6)
text(18, 1.7, "18, 1.5", cex = 0.6)
text(28, 6, "28, 5.6", cex = 0.6)

legend(0,10, legend = c('m0', 'm1_mod', 'm2'), col = c('black', 'blue', 'orange'), lty = c(2, 1, 2), box.lwd = '', cex = 0.6)
title('this work')

plot(temp_C_old, qhat_old(temp_C_old),xlim = c(0,60), ylim=c(0, 10), type ="l", xlab = '',
     ylab = expression(q['max,opt']~(g[COD-S]~g[COD-B]^'-1'~d^'-1')), col = "blue", lwd = 1,lty = 1, cex.axis = 1, cex.lab = 1)
lines(temp_C2, qhat2(temp_C2), col = "orange", lwd = 1, lty = 2)
lines(temp_C3, qhat3(temp_C3), col = "palegreen4", lwd = 1,lty = 1)
lines(temp_C4, qhat4(temp_C4), col = "purple", lwd = 1, lty = 2)
lines(temp_C5, qhat5(temp_C5), col = "red", lwd = 1, lty = 1)

text(18, 4.0, "18, 3.6", cex = 0.6)
text(28, 6.0, "28, 5.6", cex = 0.6)
text(36, 7.6, "36, 7.2", cex = 0.6)
text(43.75, 8.4, "44, 8", cex = 0.6)
text(55, 8.4, "55, 8", cex = 0.6)

legend(0,10, legend = c('m1','m2','m3','m4','m5'), col = c('blue', 'orange','palegreen4','purple','red'), lty = c(1, 2, 1,2,1), box.lwd = '', cex = 0.6)
title('original')

dev.off()

pdf('../figures/appendix_qhat_dynamics.pdf',  width = 16/2.54, height = 8/2.54)

par(mfrow=c(1,2))
par(mar=c(4,5,1,1))
par(cex.axis=1, cex.lab=1, cex.main=0.6, cex.sub=0.6)

plot(temp_C, qhat(temp_C),xlim = c(0,60), ylim=c(0, 10), type ="l", xlab = '',
     ylab = expression(q['max,opt']~(g[COD-S]~g[COD-B]^'-1'~d^'-1')), col = "blue", lwd = 1, lty = 1)
lines(temp_C0, qhat0(temp_C0), col = "black", lwd = 1, lty = 2)
lines(temp_C2, qhat2(temp_C2), col = "orange", lwd = 1, lty = 2)
text(18, 4.0, "18, 3.6", cex = 0.6)
text(18, 1.7, "18, 1.5", cex = 0.6)
text(28, 6, "28, 5.6", cex = 0.6)

legend(0,10, legend = c('m0', 'm1_mod', 'm2'), col = c('black', 'blue', 'orange'), lty = c(2, 1, 2), box.lwd = '', cex = 0.6)
title('this work')

plot(temp_C_old, qhat_old(temp_C_old),xlim = c(0,60), ylim=c(0, 10), type ="l", xlab = '',
     ylab = expression(q['max,opt']~(g[COD-S]~g[COD-B]^'-1'~d^'-1')), col = "blue", lwd = 1,lty = 1, cex.axis = 1, cex.lab = 1)
lines(temp_C2, qhat2(temp_C2), col = "orange", lwd = 1, lty = 2)
lines(temp_C3, qhat3(temp_C3), col = "palegreen4", lwd = 1,lty = 1)
lines(temp_C4, qhat4(temp_C4), col = "purple", lwd = 1, lty = 2)
lines(temp_C5, qhat5(temp_C5), col = "red", lwd = 1, lty = 1)

text(18, 4.0, "18, 3.6", cex = 0.6)
text(28, 6.0, "28, 5.6", cex = 0.6)
text(36, 7.6, "36, 7.2", cex = 0.6)
text(43.75, 8.4, "44, 8", cex = 0.6)
text(55, 8.4, "55, 8", cex = 0.6)

legend(0,10, legend = c('m1','m2','m3','m4','m5'), col = c('blue', 'orange','palegreen4','purple','red'), lty = c(1, 2, 1,2,1), box.lwd = '', cex = 0.6)
title('original')

dev.off()


