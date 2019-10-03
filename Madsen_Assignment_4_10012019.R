## Crowley et al. 2019
## 09/05/2019

## First define the BetaVec function
BetaVec <- function(Fx,a,b,res){
  p = rep(0,res) ## p vector
  q = rep(0,res) ## psi vector
  x = seq(from = 0.5/res, to = (res-0.5)/res, by = 1/res)
  bet = (x^(a-1))*((1-x)^(b-1))
  summer = 0
  for(i in 1:res){
    summer = summer +bet[i]/res
  }
  beta = bet/summer
  better = 0
  batter = 0
  for(i in res:1){
    better = better + beta[i]/res
    batter = batter + (beta[i]/res)*(i-0.5)/res
    p[i] = better
    q[i] = (batter/p[i])*Fx
  }
  return(list(beta, p, q)) 
}

## Now define the FDraw function
FDraw <- function(Cumbeta, Fx, res){
  frac = runif(n = 1, min = 0, max = 1)
  i = 1
  while(Cumbeta[i]<frac && i < res){
    i = i + 1
  }
  if(i == 1){
    Cmin = 0
  }else{
    Cmin = Cumbeta[i-1]
  }
  ffrac = i/res - (1/res)*(Cumbeta[i]-frac)/(Cumbeta[i]-Cmin)
  FF = ffrac*Fx
  return(FF)
}


## Define variables
Fx = 100
a = 4
b = 4
n = 30
h = 0.5
s = 0.9
del = 3
c = 1

FxH = 75
aH  = a
bH = b
nH = n
hH  = 0.25
sH = s
delH = del
cH = c

res = 1000
m = 1000
seePre = 1
seeMix = 1
seePost = 1

ETS = matrix(0, nrow = 1, ncol = 3)
ETS

## Calculate functions for HIREC
stuff <- BetaVec(Fx,a,b,res)
beta = stuff[[1]]
p = stuff[[2]]
q = stuff[[3]]

I = rep(0,n)
f = rep(0,n)
P = rep(0,n)

for(i in n:1){
  prod = 1
  if(i <= (n-1)){
    for(j in i:(n-1)){
      if(j > i){
        prod = prod*(1-h*P[j])
      }
      step = (s^(j-i+1))*h*P[j+1]*(I[j+1]-(((j+1)^c) - (i^c))*del)*prod
      print(step)
      if(step < 0){
        step = 0
      }
      f[i] = f[i] + step
    }
  }
  frac = f[i]/Fx
  if(frac > 1){
    frac = 1
  }
  ii = ceiling(res*frac)
  if(ii == 0){
    ii = 1
  }
  jj = res*frac - floor(res*frac)
  P[i] = p[ii] - (beta[ii]/res)*jj
  I[i] = q[ii] - (beta[ii]/res)*((ii-0.5)/res)*Fx*jj
}

prod = 1
F00 = 0
for(j in 0:(n-1)){
  if(j > 0){
    prod = prod*(1-h*P[j])
  }
  step = s^(j+1)*h*P[j+1]*(I[j+1]-((j+1)^c)*del)*prod
  if(step < 0){
    step = 0
  }
  F00 = F00+step
}
ETS[1] = 1
prod = 1
for(kk in 2:nH){
  prod = prod*(1-h*P[kk-1])
  ETS[1] = ETS[1] + s^kk*prod
}

Sim = matrix(0, nrow = m, ncol = 2)
Cumbeta = seq(1,res)
Cumbeta[1] = beta[1]/res
for(i in 2:res){
  Cumbeta[i] = Cumbeta[i-1] + beta[i]/res
}
for(mm in 1:m){
  flag = 1
  i = 1
  while (flag&&(i<=n)) {
    mort = runif(n = 1, min = 0, max = 1)
    if(mort > s){
      Sim[mm,1] = 0
      Sim[mm,2] = i
      flag = 0
    }else{
      patch = runif(n = 1, min = 0, max = 1)
      if(patch < h){
        FF = FDraw(Cumbeta, Fx, res)
        if(FF >= f[i]){
          Sim[mm, 1] = FF-(i^c)*del
          Sim[mm, 2] = i
          flag = 0
        }
      }
      
    }
    i = i + 1
  }
  if(i > n && flag != 0){
    Sim[mm, 1] = 0
    Sim[mm, 2] = n
  } 
}

FitmeanAndStepmean = colMeans(Sim)
EmFitMean = FitmeanAndStepmean[1]
EmStopMean = FitmeanAndStepmean[2]

#####
## Calculate functions for the adapted post-HIREC world
stuffH <- BetaVec(FxH,aH,bH,res)
betaH <- stuffH[[1]]
pH <- stuffH[[2]]
qH <- stuffH[[3]]

FH = seq(1,nH)
PH = seq(1,nH)
IH = seq(1,nH)

for(i in nH:1){
  prod = 1
  for(j in i:nH-1){
    if(j > i){
      prod = prod*(1 - hH*PH[j])
      }
      step = (sH^(j-i+1))*hH*PH[j+1]*(IH[j+1]-(((j+1)^cH) - (i^cH))*delH)
      if(step < 0){
        step = 0
      }
      FH[i] = FH[i] + step
    }
    frac = FH[i]/FxH
    if(frac > 1){
      frac = 1
    }
    ii = ceiling(res*frac)
    if(ii == 0){
      ii = 1
  }
  jj = res*frac - floor(res*frac)
  PH[i] = pH[ii] - (betaH[ii]/res)*jj
  IH[i] = qH[ii] - (betaH[ii]/res)*((ii-0.5)/res)*FxH*jj
}

prod = 1
F00H = 0
for(j in 0:(nH-1)){
  if(j < 0){
    prod = prod*(1-hH*PH[j])
  }
  step = (sH^(j+1))*hH*PH[j+1]*(IH[j+1]-((j+1)^cH)*delH)*prod
  if(step < 0){
    step = 0
  }
  F00H = F00H + step
}

ETS[3] = 1
prod = 1
for(kk in 2:nH){
  prod = prod*(1-hH*PH[kk-1])
  ETS[3] = ETS[3] + sH^kk*prod
}

## Pre-HIREC and Post-HIREC fitness distributions
bet = rep(0,res)
betH = rep(0,res)
beta = rep(0,res)
betaH = rep(0,res)
x = seq(from = 0.5/res, to = (res-0.5)/res, by = 1/res)
xx = x*Fx
yy = x*FxH
bet = (x^(a-1))*((1-x)^(b-1))
betH = (x^(aH-1))*((1-x)^(bH-1))
sum = 0
sumH = 0
for(i in 1:res){
  sum = sum + bet[i]/res
  sumH = sumH + betH[i]/res
}
bet = bet/sum
betaH = betH/sumH
Feta = beta/Fx
FetaH = betaH/FxH
maxF = max(Feta)
maxFH = max(FetaH)
if(maxF > maxFH){
  maxx = maxF
}else{
  maxx = maxFH
}

plot(xx, Feta, col = "black", xlim = c(xmin = 0, xmax = Fx), ylim = c(ymin = 0, ymax = maxx), xlab = '', ylab = '')
par(new = TRUE)
plot(yy, FetaH, col = "red", xlim = c(xmin = 0, xmax = Fx), ylim = c(ymin = 0, ymax = maxx), xlab = 'Patch Fitness', ylab = 'Frequency density of patch fitnesses')

## Calculate functions with pre-HIREC thresholds but post-HIREC parameters 
PX = seq(1,nH)
IX = seq(1,nH)
for(i in 1:nH){
  frac = f[i]/FxH
  if(frac > 1){
    frac = 1
  }
  ii = ceiling(res*frac)
  if(ii == 0){
    ii = 1
  }
  jj = res*frac-floor(res*frac)
  PX[i] = pH[ii]-(betaH[ii]/res)*jj
  IX[i] = qH[ii]-(betaH[ii]/res)*((ii-0.5)/res)*FxH*jj
}
prod = 1
F00X = 0
for(j in 0:(nH-1)){
  if(j > 0){
    prod = prod*(1-hH*PX[j])
  }
  step = sH^(j+1)*hH*PX[j+1]*(IX[j+1]-((j+1)^cH)*delH)*prod
  if(step < 0){
    step = 0
  }
  F00X = F00X + step
}

ETS[2] = 1
prod = 1
for(kk in 2:nH){
  prod = prod*(1-hH*PX[kk-1])
  ETS[2] = ETS[2] + sH^kk*prod
}
## simulation pre/post 
SimX = matrix(0, m, 2)
Cumbeta = seq(1, res)
Cumbeta[1] = betaH[1]/res
for(i in 2:res){
  Cumbeta[i] = Cumbeta[i-1] + betaH[i]/res
}
for(mm in 1:m){
  flag = 1
  i =1
  while (flag&&(i<=nH)) {
    mort = runif(n = 1, min = 0, max = 1)
    if(mort > sH){
      SimX[mm,1] = 0
      SimX[mm,2] = i
      flag = 0
    }else{
      patch = runif(n = 1, min = 0, max =1)
      if(patch < hH){
        FF = FDraw(Cumbeta, FxH, res)
        if(FF >= f[i]){
          SimX[mm,1] = FF - (i^cH)*delH
          SimX[mm,2] = i
          flag = 0
        }
      }
    }
    i = i + 1
  }
  if((i > nH)&&flag){
    SimX[mm,1] = 0
    SimX[mm,2] = nH
  }
}
FitmeanAndStepmeanX = colMeans(SimX)
EmFitMeanX = FitmeanAndStepmeanX[1]
EmStopMeanX = FitmeanAndStepmeanX[2]

## SImulate with post-HIREC parameters and thresholds
SimH = matrix(0, m, 2)
for(mm in 1:m){
  flag = 1
  i = 1
  while (flag&&(i<=nH)) {
    mort = runif(n = 1, min = 0, max = 1)
    if(mort > sH){
      SimH[mm,1] = 0
      SimH[mm,2] = i
      flag = 0
    }else{
      patch = runif(n = 1, min = 0, max = 1)
      if(patch < hH){
        FFH = FDraw(Cumbeta, FxH, res)
        if(FFH >= FH[i]){
          SimH[mm,1] = FFH-(i^cH)*delH
          SimH[mm,2] = i
          flag = 0
        }
      }
    }
    i = i +1
  }
  if((i > nH)&&flag){
    SimH[mm,1] = 0
    SimH[mm,2] = n
  }
}
FitmeanAndStepmeanH = colMeans(SimH)
EmFitMeanH = FitmeanAndStepmeanH[1]
EmStopMeanH = FitmeanAndStepmeanH[2]

if(nH > n){
  nmx = nH
}else{
  nmx = n
}
if(FxH > Fx){
  Fxmx = FxH
}else{
  Fxmx = Fx
}
## Output plot
xx = 1:n
xxx = 1:nH
plot(xx, f, col = "black", xlim = c(xmin = 0, xmax = 30), ylim = c(ymin = 0, ymax = 100), pch = 4, xlab = '', ylab = '')
par(new = TRUE)
plot(xxx, FH, xlim = c(xmin = 0, xmax = 30), ylim = c(ymin = 0, ymax = 100), col = "blue", pch = 4, xlab = '', ylab = '')
abline(h = F00, col = "black", pch = 19)
abline(h = F00X, col = "red", pch = 19)
abline(h = F00H, col = "blue", pch = 19)
par(new = TRUE)
plot(Sim[,2], Sim[,1], xlim = c(xmin = 0, xmax = 30), ylim = c(ymin = 0, ymax = 100), col = "black", pch = 20, xlab = 'Time Step', ylab = 'Fitness')
par(new = TRUE)
plot(SimX[,2], SimX[,1], xlim = c(xmin = 0, xmax = 30), ylim = c(ymin = 0, ymax = 100), col = "red", pch = 20, xlab = '', ylab = '')
par(new = TRUE)
plot(SimH[,2], SimH[,1], xlim = c(xmin = 0, xmax = 30), ylim = c(ymin = 0, ymax = 100), col = "blue", pch = 20, xlab = '', ylab = '')
