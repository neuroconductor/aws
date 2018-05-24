TGV_denoising <- function(alpha,beta,iter,datanoisy){
# translated from original matlab code of K. Papafitsoros
# Noisy data
f <- datanoisy
df <- dim(f)
n <- df[1]
m <- df[2]
# TGV parameters (TGV(u)=min_{w} a1|| nabla u-w ||_{1}+a0|| Ew ||_{1})
a0 <- beta
a1 <- alpha
N <- iter
# Initializations
ubar <- uold <- f
pbar1 <- p1_old <- pbar2 <- p2_old  <-
v1 <- v2 <- w11 <- w12 <- w22 <- matrix(0,n,m)
# Setting the parameters sigma and tau of Chambolle-Pock method
L2 <- 8

for (k in 1:N){
#    Acceleration for the parameters sigma and tau
    tau <- 1/(k+1)
    sigma <- (1/tau)/L2

    v1 <- v1+sigma*(Dx_plus(ubar)-pbar1)
    v2 <- v2+sigma*(Dy_plus(ubar)-pbar2)
    z <- pmax(1,sqrt(v1^2+v2^2)/a1)
    v1 <- v1/z
    v2 <- v2/z

    w11 <- w11+sigma*(Dx_minus(pbar1))
    w22 <- w22+sigma*(Dy_minus(pbar2))
    w12 <- w12+sigma*(0.5*(Dy_minus(pbar1)+Dx_minus(pbar2)))
    z <- pmax(1,sqrt(w11^2+w22^2+w12^2)/a0)
    w11 <- w11/z
    w12 <- w12/z
    w22 <- w22/z

    unew <- (1/(1+tau))*(uold+tau*(Dx_minus(v1)+Dy_minus(v2)+f))

    p1_new <- p1_old+tau*(v1+Dx_plus(w11)+Dy_plus(w12))
    p2_new <- p2_old+tau*(v2+Dy_plus(w22)+Dx_plus(w12))

    ubar <- 2*unew-uold

    pbar1 <- 2*p1_new-p1_old
    pbar2 <- 2*p2_new-p2_old
    uold <- unew

    p1_old <- p1_new
    p2_old <- p2_new
}
cat(iter,"Chambolle-Pock iterations completed\n")
unew
}

TV_denoising <- function(alpha,iter,datanoisy){
# Noisy data
f <- datanoisy
df <- dim(f)
n <- df[1]
m <- df[2]

# TV parameter: We are solving min_{u} 0.5|| u-f ||_{2}^{2}+alpha TV(u)
# or equivalently              min_{u} 0.5lambda|| u-f ||_{2}^{2}+TV(u)
# with                         lambda=1/alpha
lambda <- 1/alpha
N <- iter

# Initializations
ubar <- u_old <- f
p1 <- p2 <- matrix(0,n,m)
# Setting the parameters sigma and tau of Chambolle-Pock method
L2 <- 8

for (k in 1:N) {
#  Acceleration for the parameters sigma and tau
    tau <- 1/(k+1)
    sigma <- (1/tau)/L2
    ltau <- lambda*tau

    p1 <- p1+sigma*(Dx_plus(ubar))
    p2 <- p2+sigma*(Dy_plus(ubar))
    z <- pmax(1,sqrt(p1^2+p2^2))
    p1 <- p1/z
    p2 <- p2/z

    u_new <- (1/(1+ltau))*(u_old+tau*
                (Dx_minus(p1)+Dy_minus(p2))+ltau*f)

    ubar <- 2*u_new-u_old
    u_old <- u_new
  }
  cat(iter,"Chambolle-Pock iterations completed\n")
  u_new
}

TGV_denoising_colour <- function(alpha,beta,iter,datanoisy){
# translated from original matlab code of K. Papafitsoros
# Noisy data
f <- datanoisy
f_r <- f[,,1]
f_g <- f[,,2]
f_b <- f[,,3]

df <- dim(f)
n <- df[1]
m <- df[2]
# TGV parameters (TGV(u)=min_{w} a1|| nabla u-w ||_{1}+a0|| Ew ||_{1})
a0 <- beta
a1 <- alpha
N <- iter
# Initializations
ubar_r <- uold_r <- f_r
ubar_g <- uold_g <- f_g
ubar_b <- uold_b <- f_b

pbar1_r <- p1_old_r <- pbar2_r <- p2_old_r <-
v1_r <- v2_r <- w11_r <- w12_r <- w22_r <- matrix(0,n,m)

pbar1_g <- p1_old_g <- pbar2_g <- p2_old_g <-
v1_g <- v2_g <- w11_g <- w12_g <- w22_g <- matrix(0,n,m)

pbar1_b <- p1_old_b <- pbar2_b <- p2_old_b <-
v1_b <- v2_b <- w11_b <- w12_b <- w22_b <- matrix(0,n,m)

# Setting the parameters sigma and tau of Chambolle-Pock method
L2 <- 8

for (k in 1:N){
#    Acceleration for the parameters sigma and tau
    tau <- 1/(k+1)
    sigma <- (1/tau)/L2

    v1_r <- v1_r+sigma*(Dx_plus(ubar_r)-pbar1_r)
    v2_r <- v2_r+sigma*(Dy_plus(ubar_r)-pbar2_r)

    v1_g <- v1_g+sigma*(Dx_plus(ubar_g)-pbar1_g)
    v2_g <- v2_g+sigma*(Dy_plus(ubar_g)-pbar2_g)

    v1_b <- v1_b+sigma*(Dx_plus(ubar_b)-pbar1_b)
    v2_b <- v2_b+sigma*(Dy_plus(ubar_b)-pbar2_b)

    z <- pmax(1,sqrt(v1_r^2+v2_r^2+v1_g^2+v2_g^2+v1_b^2+v2_b^2)/a1)
    v1_r <- v1_r/z
    v2_r <- v2_r/z
    v1_g <- v1_g/z
    v2_g <- v2_g/z
    v1_b <- v1_b/z
    v2_b <- v2_b/z

    w11_r <- w11_r+sigma*(Dx_minus(pbar1_r))
    w22_r <- w22_r+sigma*(Dy_minus(pbar2_r))
    w12_r <- w12_r+sigma*(0.5*(Dy_minus(pbar1_r)+Dx_minus(pbar2_r)))

    w11_g <- w11_g+sigma*(Dx_minus(pbar1_g))
    w22_g <- w22_g+sigma*(Dy_minus(pbar2_g))
    w12_g <- w12_g+sigma*(0.5*(Dy_minus(pbar1_g)+Dx_minus(pbar2_g)))

    w11_b <- w11_b+sigma*(Dx_minus(pbar1_b))
    w22_b <- w22_b+sigma*(Dy_minus(pbar2_b))
    w12_b <- w12_b+sigma*(0.5*(Dy_minus(pbar1_b)+Dx_minus(pbar2_b)))

    z <- pmax(1,sqrt(w11_r^2+w12_r^2+w22_r^2+
                     w11_g^2+w12_g^2+w22_g^2+
                     w11_b^2+w12_b^2+w22_b^2)/a0)
    w11_r <- w11_r/z
    w12_r <- w12_r/z
    w22_r <- w22_r/z

    w11_g <- w11_g/z
    w12_g <- w12_g/z
    w22_g <- w22_g/z

    w11_b <- w11_b/z
    w12_b <- w12_b/z
    w22_b <- w22_b/z

    unew_r <- (1/(1+tau))*(uold_r+tau*(Dx_minus(v1_r)+Dy_minus(v2_r)+f_r))
    p1_new_r <- p1_old_r+tau*(v1_r+Dx_plus(w11_r)+Dy_plus(w12_r))
    p2_new_r <- p2_old_r+tau*(v2_r+Dy_plus(w22_r)+Dx_plus(w12_r))

    unew_g <- (1/(1+tau))*(uold_g+tau*(Dx_minus(v1_g)+Dy_minus(v2_g)+f_g))
    p1_new_g <- p1_old_g+tau*(v1_g+Dx_plus(w11_g)+Dy_plus(w12_g))
    p2_new_g <- p2_old_g+tau*(v2_g+Dy_plus(w22_g)+Dx_plus(w12_g))

    unew_b <- (1/(1+tau))*(uold_b+tau*(Dx_minus(v1_b)+Dy_minus(v2_b)+f_b))
    p1_new_b <- p1_old_b+tau*(v1_b+Dx_plus(w11_b)+Dy_plus(w12_b))
    p2_new_b <- p2_old_b+tau*(v2_b+Dy_plus(w22_b)+Dx_plus(w12_b))

    ubar_r <- 2*unew_r-uold_r
    pbar1_r <- 2*p1_new_r-p1_old_r
    pbar2_r <- 2*p2_new_r-p2_old_r

    ubar_g <- 2*unew_g-uold_g
    pbar1_g <- 2*p1_new_g-p1_old_g
    pbar2_g <- 2*p2_new_g-p2_old_g

    ubar_b <- 2*unew_b-uold_b
    pbar1_b <- 2*p1_new_b-p1_old_b
    pbar2_b <- 2*p2_new_b-p2_old_b

    relchange <- abs(uold_r-unew_r)+abs(uold_g-unew_g)+abs(uold_b-unew_b)
    maxrelchange <- max(relchange)
    meanrelchange <- mean(relchange)

    uold_r <- unew_r
    p1_old_r <- p1_new_r
    p2_old_r <- p2_new_r

    uold_g <- unew_g
    p1_old_g <- p1_new_g
    p2_old_g <- p2_new_g

    uold_b <- unew_b
    p1_old_b <- p1_new_b
    p2_old_b <- p2_new_b

    cat(k," out of ",iter," iterations completed, changed (mean,max)",meanrelchange,maxrelchange,"\n")
}
cat(iter,"Chambolle-Pock iterations completed\n")
f[,,1] <- unew_r
f[,,2] <- unew_g
f[,,3] <- unew_b
f
}

TV_denoising_colour <- function(alpha,iter,datanoisy){
# Noisy data
f <- datanoisy
f_r <- f[,,1]
f_g <- f[,,2]
f_b <- f[,,3]
df <- dim(f)
n <- df[1]
m <- df[2]

# TV parameter: We are solving min_{u} 0.5|| u-f ||_{2}^{2}+alpha TV(u)
# or equivalently              min_{u} 0.5lambda|| u-f ||_{2}^{2}+TV(u)
# with                         lambda=1/alpha
lambda <- 1/alpha
N <- iter

# Initializations
ubar_r <- u_old_r <- f_r
p1_r <- p2_r <- matrix(0,n,m)

ubar_g <- u_old_g <- f_g
p1_g <- p2_g <- matrix(0,n,m)

ubar_b <- u_old_b <- f_b
p1_b <- p2_b <- matrix(0,n,m)

# Setting the parameters sigma and tau of Chambolle-Pock method
L2 <- 8

for (k in 1:N) {
#  Acceleration for the parameters sigma and tau
    tau <- 1/(k+1)
    sigma <- (1/tau)/L2;

    p1_r <- p1_r+sigma*(Dx_plus(ubar_r))
    p2_r <- p2_r+sigma*(Dy_plus(ubar_r))

    p1_g <- p1_g+sigma*(Dx_plus(ubar_g))
    p2_g <- p2_g+sigma*(Dy_plus(ubar_g))

    p1_b <- p1_b+sigma*(Dx_plus(ubar_b))
    p2_b <- p2_b+sigma*(Dy_plus(ubar_b))

    z <- pmax(1,sqrt(p1_r^2+p2_r^2+p1_g^2+p2_g^2+p1_b^2+p2_b^2))
    p1_r <- p1_r/z
    p2_r <- p2_r/z
    p1_g <- p1_g/z
    p2_g <- p2_g/z
    p1_b <- p1_b/z
    p2_b <- p2_b/z

    u_new_r <- (1/(1+lambda*tau))*(u_old_r+tau*
                  (Dx_minus(p1_r)+Dy_minus(p2_r))+lambda*tau*f_r)
    u_new_g <- (1/(1+lambda*tau))*(u_old_g+tau*
                  (Dx_minus(p1_g)+Dy_minus(p2_g))+lambda*tau*f_g)
    u_new_b <- (1/(1+lambda*tau))*(u_old_b+tau*
                  (Dx_minus(p1_b)+Dy_minus(p2_b))+lambda*tau*f_b)

    ubar_r <- 2*u_new_r-u_old_r
    ubar_g <- 2*u_new_g-u_old_g
    ubar_b <- 2*u_new_b-u_old_b

    relchange <- abs(u_old_r-u_new_r)+abs(u_old_g-u_new_g)+abs(u_old_b-u_new_b)
    maxrelchange <- max(relchange)
    meanrelchange <- mean(relchange)

    u_old_r <- u_new_r
    u_old_g <- u_new_g
    u_old_b <- u_new_b
    cat(k," out of ",iter," iterations completed, changed (mean,max)",meanrelchange,maxrelchange,"\n")
}
  cat(iter,"Chambolle-Pock iterations completed\n")
  f[,,1] <- u_new_r
  f[,,2] <- u_new_g
  f[,,3] <- u_new_b
  f
}

P_a0 <- function(w11,w22,w12,a0){
normw_over_a0 <- sqrt(w11^2+w22^2+2*w12^2)/a0
denom= pmax(1,normw_over_a0)
list(w11=w11/denom, w22=w22/denom, w12=w12/denom)
}

P_a1 <- function(v1,v2,a1){
normv_over_a1 <- sqrt(v1^2+v2^2)
denom=pmax(1,normv_over_a1/a1)
list(v1=v1/denom, v2=v2/denom)
}
dP_a0 <- function(w11,w22,w12,a0){
normw_over_a0 <- sqrt(w11^2+w22^2+2*w12^2)/a0
pmax(1,normw_over_a0)
}

dP_a1 <- function(v1,v2,a1){
normv_over_a1 <- sqrt(v1^2+v2^2)
pmax(1,normv_over_a1/a1)
}

P_a0_colour <- function(w11_r,w11_g,w11_b,w22_r,w22_g,w22_b,w12_r,w12_g,w12_b,a0){
normw_over_a0 <- (1/a0)*sqrt(w11_r^2+w22_r^2+2*w12_r^2+
                             w11_g^2+w22_g^2+2*w12_g^2+
                             w11_b^2+w22_b^2+2*w12_b^2)
denom= pmax(1,normw_over_a0)
list(w11_r=w11_r/denom, w22_r=w22_r/denom, w12_r=w12_r/denom,
     w11_g=w11_g/denom, w22_g=w22_g/denom, w12_g=w12_g/denom,
     w11_b=w11_b/denom, w22_b=w22_b/denom, w12_b=w12_b/denom)
}

P_a1_colour <- function(v1_r,v2_r,v1_g,v2_g,v1_b,v2_b,a1){
normv_over_a1=sqrt(v1_r^2+v2_r^2+v1_g^2+v2_g^2+v1_b^2+v2_b^2)
denom=pmax(1,normv_over_a1/a1)
list(v1_r=v1_r/denom, v2_r=v2_r/denom,v1_g=v1_g/denom,
     v2_g=v2_g/denom,v1_b=v1_b/denom, v2_b=v2_b/denom)
}

dP_a0_colour <- function(w11_r,w11_g,w11_b,w22_r,w22_g,w22_b,w12_r,w12_g,w12_b,a0){
normw_over_a0 <- sqrt(w11_r^2+w22_r^2+2*w12_r^2+
                      w11_g^2+w22_g^2+2*w12_g^2+
                      w11_b^2+w22_b^2+2*w12_b^2)/a0
pmax(1,normw_over_a0)
}

dP_a1_colour <- function(v1_r,v2_r,v1_g,v2_g,v1_b,v2_b,a1){
normv_over_a1 <- sqrt(v1_r^2+v2_r^2+v1_g^2+v2_g^2+v1_b^2+v2_b^2)
pmax(1,normv_over_a1/a1)
}

Dx_minus <- function(u){
  du <- dim(u)
  n <- du[1]
  m <- du[2]
v <- u
v[2:(n-1),] <- u[2:(n-1),] - u[1:(n-2),]
v[n,] <- -u[n-1,]
v
}
Dx_plus <- function(u){
du <- dim(u)
n <- du[1]
m <- du[2]
v <- u
v[1:(n-1),] <- u[2:n,] - u[1:(n-1),]
v[n,] <- 0
v
}
Dy_minus <- function(u){
  du <- dim(u)
  n <- du[1]
  m <- du[2]
v <- u
v[,2:(m-1)] <- u[,2:(m-1)] - u[,1:(m-2)]
v[,m] <- -u[,m-1]
v
}
Dy_plus <- function(u){
  du <- dim(u)
  n <- du[1]
  m <- du[2]
v <- u
v[,1:(m-1)] <- u[,2:m] - u[,1:(m-1)]
v[,m] <- 0
v
}
