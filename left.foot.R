library('RMoCap')

#set your path
path <- "e:\\Publikacje\\Danuta kata\\gotowy_kod\\" 

shodan.marcin <- read.mocap(paste(path, "data\\pinian_shodan.bvh", sep = ""))
shodan.danuta <- read.mocap(paste(path, "data\\1-pinian_shodan.bvh", sep = ""))

#uncomment if you want to use other pairs
#shodan.marcin <- read.mocap(paste(path, "data\\pinian_nidan.bvh", sep = ""))
#shodan.danuta <- read.mocap(paste(path, "data\\2-pinian_nidan.bvh", sep = ""))

#shodan.marcin <- read.mocap(paste(path, "data\\fukyugata_ni.bvh", sep = ""))
#shodan.danuta <- read.mocap(paste(path, "data\\0.2-fukyugata_ni.bvh")#kata3.bvh", sep = ""))

#shodan.marcin <- read.mocap(paste(path, "data\\fukyugata_ichi.bvh", sep = ""))
#shodan.danuta <- read.mocap(paste(path, "data\\0.1-fukyugata_ichi.bvh")#kata3.bvh", sep = ""))


input <- shodan.marcin
ref <- shodan.danuta

pca.input.x <- input$data.frame$LeftFoot.Dx
pca.input.z <- input$data.frame$LeftFoot.Dz

pca.ref.x <- ref$data.frame$LeftFoot.Dx
pca.ref.z <- ref$data.frame$LeftFoot.Dz

plot(x = pca.input.x, y = pca.input.z, type = "l", col = "blue", 
     xlim = c(min(pca.input.x, pca.ref.x),
              max(pca.input.x, pca.ref.x)),
     ylim = c(min(pca.input.z, pca.ref.z),
              max(pca.input.z, pca.ref.z)), xlab="X-axis [cm]", ylab="y-axis [cm]")
lines(x = pca.ref.x, y = pca.ref.z, col = "green")

points(x = pca.input.x[1], y = pca.input.z[1], col = "blue", pch = 1, lwd = 5)
points(x = pca.ref.x[1], y = pca.ref.z[1], col = "green", pch = 5, lwd = 5)

points(x = pca.input.x[length(pca.input.x)], y = pca.input.z[length(pca.input.x)], col = "blue", pch = 1, lwd = 5)
points(x = pca.ref.x[length(pca.ref.x)], y = pca.ref.z[length(pca.ref.z)], col = "green", pch = 5, lwd = 5)


input.raw <- input$skeleton$Joints[[1]]$RawDxyz[1,]
for (a in 1:length(input$skeleton$Joints[[1]]$Rxyz[,1]))
{
  input$skeleton$Joints[[1]]$RawDxyz[a, 1] <- input$skeleton$Joints[[1]]$RawDxyz[a,1] - input.raw[1]
  input$skeleton$Joints[[1]]$RawDxyz[a, 3] <- input$skeleton$Joints[[1]]$RawDxyz[a,3] - input.raw[3]
}
input$data.frame <- hierarchical.to.direct.kinematic(input$skeleton)


ref.raw <- ref$skeleton$Joints[[1]]$RawDxyz[1,]
for (a in 1:length(ref$skeleton$Joints[[1]]$Rxyz[,1]))
{
  ref$skeleton$Joints[[1]]$RawDxyz[a, 1] <- ref$skeleton$Joints[[1]]$RawDxyz[a,1] - ref.raw[1]
  ref$skeleton$Joints[[1]]$RawDxyz[a, 3] <- ref$skeleton$Joints[[1]]$RawDxyz[a,3] - ref.raw[3]
}
ref$data.frame <- hierarchical.to.direct.kinematic(ref$skeleton)



###########################################

library('subplex')

input <- shodan.marcin
ref <- shodan.danuta

input.corrected <- aligninputandrefdata(refdata = ref$data.frame,inputdata = input$data.frame, limbname = "Hips")

input$skeleton$Joints[[1]]$RawDxyz[,1] <- input.corrected$Hips.Dx
input$skeleton$Joints[[1]]$RawDxyz[,2] <- input.corrected$Hips.Dy
input$skeleton$Joints[[1]]$RawDxyz[,3] <- input.corrected$Hips.Dz

input$data.frame <- hierarchical.to.direct.kinematic(input$skeleton)


v1 <- "LeftThigh"
v2 <-"RightThigh"

library('subplex')
v1.x <- paste(v1,".Dx", sep = "")
v1.z <- paste(v1,".Dz", sep = "")

v2.x <- paste(v2,".Dx", sep = "")
v2.z <- paste(v2,".Dz", sep = "")

library(RSpincalc)

euler2quaternion <- function(xR, yR, zR)
{
  a1 <- c(zR, yR, xR) * (pi/180)
  q <- EA2Q(a1,'zyx')
  return(q)
}

#helper function
quaternion2euler <- function(q)
{
  ea <- Q2EA(q,'zyx') * (180/pi)
  return(c(ea[3], ea[2], ea[1]))
}

myEV2Q <- function(axis, angle)
{
  qx = axis[1] * sin(angle/2)
  qy = axis[2] * sin(angle/2)
  qz = axis[3] * sin(angle/2)
  qw = cos(angle/2)
  Q = c(qw,qx,qy,qz)
  Q
}

#ref <- shodan.marcin
#input <- shodan.danuta
q1 <- euler2quaternion(input$skeleton$Joints[[1]]$Rxyz[1,1], input$skeleton$Joints[[1]]$Rxyz[1,2], input$skeleton$Joints[[1]]$Rxyz[1,3])





euc.dist <- function(x1, x2) sqrt(sum((x1 - x2) ^ 2))

optimizeangle <- function(x)
{
  q2 <- myEV2Q(c(0,1,0),x)
  q3 <- q2 %Q*% q1
  xyz <- quaternion2euler(q3)
  input$skeleton$Joints[[1]]$Rxyz[1,1] <- xyz[1]
  input$skeleton$Joints[[1]]$Rxyz[1,2] <- xyz[2]
  input$skeleton$Joints[[1]]$Rxyz[1,3] <- xyz[3]

  input$data.frame <- hierarchical.to.direct.kinematic(input$skeleton)
  
  vvv1 <- c(input$data.frame[1,v1.x],0, input$data.frame[1,v1.z])
  - c(input$data.frame[1,v2.x],0, input$data.frame[1,v2.z])
  
  
  vvv2 <- c(ref$data.frame[1,v1.x],0, ref$data.frame[1,v1.z])
  - c(ref$data.frame[1,v2.x],0, ref$data.frame[1,v2.z])

  message(euc.dist(vvv1, vvv2))
  return (euc.dist(vvv1, vvv2))
}

plot(ref, append = FALSE, my.color = "green", frame = 1, alpha = 1, spheres = TRUE)
plot(input, append = TRUE, my.color = "blue", frame = 1, alpha = 1, spheres = TRUE)

set.seed(1)

response <- subplex(par=c(0),fn=optimizeangle)

for (a in 1:length(input$skeleton$Joints[[1]]$Rxyz[,1]))
{
  q1 <- euler2quaternion(input$skeleton$Joints[[1]]$Rxyz[a,1], input$skeleton$Joints[[1]]$Rxyz[a,2], input$skeleton$Joints[[1]]$Rxyz[a,3])
  q2 <- myEV2Q(c(0,1,0), response$par)
  q3 <- q2 %Q*% q1
  xyz <- quaternion2euler(q3)
  
  input$skeleton$Joints[[1]]$Rxyz[a,1] <- xyz[1]
  input$skeleton$Joints[[1]]$Rxyz[a,2] <- xyz[2]
  input$skeleton$Joints[[1]]$Rxyz[a,3] <- xyz[3]
}
input$data.frame <- hierarchical.to.direct.kinematic(input$skeleton)

input$data.frame <- hierarchical.to.direct.kinematic(input$skeleton)
write.bvh(input, paste(path, "temp\\corrected_help.bvh", sep = ""))
input <- read.mocap(paste(path, "temp\\corrected_help.bvh", sep = ""))

ref$data.frame <- hierarchical.to.direct.kinematic(ref$skeleton)
write.bvh(ref, paste(path, "temp\\corrected_help.bvh", sep = ""))
ref <- read.mocap(paste(path, "temp\\corrected_help.bvh", sep = ""))


input$data.frame <- calculate.kinematic(input$data.frame, show.plot = "FALSE", plot.title = "Heian Shodan")

plot(ref, append = FALSE, my.color = "green", frame = 1, alpha = 1, spheres = TRUE)
plot(input, append = TRUE, my.color = "blue", frame = 1, alpha = 1, spheres = TRUE)

##########################################

input.dx <- input$data.frame$LeftFoot.Dx
input.dz <- input$data.frame$LeftFoot.Dz

ref.dx <- ref$data.frame$LeftFoot.Dx
ref.dz <- ref$data.frame$LeftFoot.Dz


plot(x = input.dx, y = input.dz, type = "l", col = "blue", 
     xlim = c(min(input.dx, ref.dx),
              max(input.dx, ref.dx)),
     ylim = c(min(input.dz, ref.dz),
              max(input.dz, ref.dz)), xlab="X-axis [cm]", ylab="y-axis [cm]")
lines(x = ref.dx, y = ref.dz, col = "green")

points(x = input.dx[1], y = input.dz[1], col = "blue", pch = 1, lwd = 5)
points(x = ref.dx[1], y = ref.dz[1], col = "green", pch = 5, lwd = 5)

points(x = input.dx[length(input.dx)], y = input.dz[length(input.dx)], col = "blue", pch = 1, lwd = 5)
points(x = ref.dx[length(ref.dx)], y = ref.dz[length(ref.dx)], col = "green", pch = 5, lwd = 5)


##########################################

left.right <- function(dd)
{
  #dd <- input$data.frame
  LeftFoot = "LeftFoot"
  RightFoot = "RightFoot"
  window_size <- 100 / length(dd[,paste(RightFoot,".ax", sep ="")])
  library(smoother)
  zzR <- smth(sqrt(dd[,paste(RightFoot,".ax", sep ="")] ^ 2 + dd[,paste(RightFoot,".ay", sep ="")] ^ 2 + dd[,paste(RightFoot,".az", sep ="")] ^ 2),window = window_size,method = "gaussian") #SMOOTHING
  zzL <- smth(sqrt(dd[,paste(LeftFoot,".ax", sep ="")] ^ 2 + dd[,paste(LeftFoot,".ay", sep ="")] ^ 2 + dd[,paste(LeftFoot,".az", sep ="")] ^ 2),window = window_size,method = "gaussian") #SMOOTHING
  
  zzR[is.na(zzR)] <- 0
  zzL[is.na(zzL)] <- 0

  lr <- rep(0, length(zzR))
  
  zzRy <- smth(dd[,paste(RightFoot,".Dy", sep ="")],window = window_size,method = "gaussian") #SMOOTHING
  zzLy <- smth(dd[,paste(LeftFoot,".Dy", sep ="")],window = window_size,method = "gaussian") #SMOOTHING
  
  zzRy[is.na(zzRy)] <- 0
  zzLy[is.na(zzLy)] <- 0
  
  dyEps = 5
  for (a in 1:(length(dd[,paste(RightFoot,".ax", sep ="")]) - 1))
  {
    if (zzR[a] > zzL[a] && zzR[a + 1] > zzL[a + 1])
    {
      lr[a] <- 1
    }
    else if (zzR[a] < zzL[a] && zzR[a + 1] < zzL[a + 1])
    {
      lr[a] <- -1
    }
    
    #if (abs(zzRy[a] - zzLy[a]) < dyEps)
    #{
    #  lr[a] <- 0
    #}
  }
  return (lr)
}

generate.lists <-  function(lr.input)
{
  prev <- lr.input[1]
  start <- 0
  end <- 0
  list.right <- list()
  list.left <- list()
  
  for (a in 2:length(lr.input))
  {
    if (prev != lr.input[a])
    {
      #if (lr.input[a] != 0)
      {
        if (start == 0)
        {
          start = a
        }
        else {
          end = a - 1
          if (prev == 1)
          {
            list.right[[length(list.right)+1]] <- c(start, end)
            start = 0
            end = 0
          }
          if (prev == -1)
          {
            list.left[[length(list.left)+1]] <- c(start, end)
            start = 0
            end = 0
          }
        }
      }
    }
    prev <- lr.input[a]
  }
  ret.list <- list()
  ret.list$list.left <- list.left
  ret.list$list.right <- list.right
  return(ret.list)
}

lr.input <- left.right(input$data.frame)

lr.input <- runmed(x = lr.input, k = 9)

ret.list.input <- generate.lists(lr.input)
list.left.input <- ret.list.input$list.left
list.right.input <- ret.list.input$list.right



lr.ref <- left.right(ref$data.frame)

lr.ref <- runmed(x = lr.ref, k = 9)

ret.list.ref <- generate.lists(lr.ref)
list.left.ref <- ret.list.ref$list.left
list.right.ref <- ret.list.ref$list.right

################

#input.dx <- input$data.frame$Hips.Dx
#input.dz <- input$data.frame$Hips.Dz

input.dx <- input$data.frame$LeftFoot.Dx
input.dz <- input$data.frame$LeftFoot.Dz

window_size <- 100 / length(input.dx)

options('smoother.tails'= TRUE)
library(smoother)
input.dx <- smth(input.dx,window = window_size,method = "gaussian") #SMOOTHING
input.dz <- smth(input.dz,window = window_size,method = "gaussian") #SMOOTHING


input.dx[is.na(input.dx)] <- 0
input.dz[is.na(input.dz)] <- 0

ref.dx <- ref$data.frame$LeftFoot.Dx
ref.dz <- ref$data.frame$LeftFoot.Dz

ref.dx <- smth(ref.dx,window = window_size,method = "gaussian") #SMOOTHING
ref.dz <- smth(ref.dz,window = window_size,method = "gaussian") #SMOOTHING

ref.dx[is.na(ref.dx)] <- 0
ref.dz[is.na(ref.dz)] <- 0

x0 <- c()
y0 <- c()

eps.arrow.length <- 5

list.right.input2 <- list()
i <- 1
for (a in 1:length(list.right.input))
{
  coord <- list.right.input[[a]]
  
  if ((input.dx[coord[1]] - input.dx[coord[2]]) ^ 2 + (input.dz[coord[1]] - input.dz[coord[2]]) ^ 2 > eps.arrow.length)
  {
    x0 <- c(x0, input.dx[coord[1]], input.dx[coord[2]])
    y0 <- c(y0, input.dz[coord[1]], input.dz[coord[2]])
    
    list.right.input2[[i]] <- coord
    i <- i +1
  }
}

list.right.input <- list.right.input2

list.left.input2 <- list()
i <- 1
for (a in 1:length(list.left.input))
{
  coord <- list.left.input[[a]]
  if ((input.dx[coord[1]] - input.dx[coord[2]]) ^ 2 + (input.dz[coord[1]] - input.dz[coord[2]]) ^ 2 > eps.arrow.length)
  {
    x0 <- c(x0, input.dx[coord[1]], input.dx[coord[2]])
    y0 <- c(y0, input.dz[coord[1]], input.dz[coord[2]])
    list.left.input2[[i]] <- coord
    i <- i +1
  }
}

list.left.input <- list.left.input2
################

list.right.ref2 <- list()
i <- 1
for (a in 1:length(list.right.ref))
{
  coord <- list.right.ref[[a]]
  
  if ((ref.dx[coord[1]] - ref.dx[coord[2]]) ^ 2 + (ref.dz[coord[1]] - ref.dz[coord[2]]) ^ 2 > eps.arrow.length)
  {
    x0 <- c(x0, ref.dx[coord[1]], ref.dx[coord[2]])
    y0 <- c(y0, ref.dz[coord[1]], ref.dz[coord[2]])
    
    list.right.ref2[[i]] <- coord
    i <- i +1
  }
}

list.right.ref <- list.right.ref2

list.left.ref2 <- list()

i <- 1
for (a in 1:length(list.left.ref))
{
  coord <- list.left.ref[[a]]
  
  if ((ref.dx[coord[1]] - ref.dx[coord[2]]) ^ 2 + (ref.dz[coord[1]] - ref.dz[coord[2]]) ^ 2 > eps.arrow.length)
  {
    x0 <- c(x0, ref.dx[coord[1]], ref.dx[coord[2]])
    y0 <- c(y0, ref.dz[coord[1]], ref.dz[coord[2]])
    
    list.left.ref2[[i]] <- coord
    i <- i +1
  }
}

list.left.ref <- list.left.ref2

######################

plot(x = input.dx, y = input.dz, type = "l", col = "blue", 
      xlim = c(min(input.dx, ref.dx),
               max(input.dx, ref.dx)),
      ylim = c(min(input.dz, ref.dz),
               max(input.dz, ref.dz)), lwd=1, xlab="X-axis [cm]", ylab="y-axis [cm]")
lines(x = ref.dx, y = ref.dz, col = "green", lwd=1)

ll1 <- list()
ml1 <- matrix(rep(0, 2 * length(list.left.input)), nrow = length(list.left.input), ncol = 2)
for (a in 1:length(list.left.input))
{
  coord <- list.left.input[[a]]

  ll1[[a]] <- c(input.dx[coord[2]] - input.dx[coord[1]], input.dz[coord[2]] - input.dz[coord[1]])
  ml1[a, ] <- c(input.dx[coord[2]] - input.dx[coord[1]], input.dz[coord[2]] - input.dz[coord[1]])
  arrows(x0 = input.dx[coord[1]], y0 = input.dz[coord[1]], x1 = input.dx[coord[2]], y1 = input.dz[coord[2]], 
         length = 0.25, angle = 30, col = 'blue', lwd=4)
}



ll2 <- list()
ml2 <- matrix(rep(0, 2 * length(list.left.ref)), nrow = length(list.left.ref), ncol = 2)
for (a in 1:length(list.left.ref))
{
  coord <- list.left.ref[[a]]
  
  ll2[[a]] <- c(ref.dx[coord[2]] - ref.dx[coord[1]], ref.dz[coord[2]] - ref.dz[coord[1]])
  ml2[a, ] <- c(ref.dx[coord[2]] - ref.dx[coord[1]], ref.dz[coord[2]] - ref.dz[coord[1]])
  arrows(x0 = ref.dx[coord[1]], y0 = ref.dz[coord[1]], x1 = ref.dx[coord[2]], y1 = ref.dz[coord[2]], 
         length = 0.25, angle = 30, col = 'green', lwd=4)
}

library(dtw)
my.dtw <- dtw(x = ml1, y = ml2, keep=TRUE)

for (a in 1:length(my.dtw$index1))
{
  coord.in <- list.left.input[[my.dtw$index1[a]]]
  coord.ref <- list.left.ref[[my.dtw$index2[a]]]

  lines(x = c((input.dx[coord.in[1]] + input.dx[coord.in[2]]) / 2, (ref.dx[coord.ref[1]] + ref.dx[coord.ref[2]]) / 2), 
             y = c((input.dz[coord.in[1]] + input.dz[coord.in[2]]) / 2, (ref.dz[coord.ref[1]] + ref.dz[coord.ref[2]]) / 2), 
        col = "red", lty = 3, lwd=2)
  
  text(x = ((input.dx[coord.in[1]] + input.dx[coord.in[2]]) / 2), 
       y = ((input.dz[coord.in[1]] + input.dz[coord.in[2]]) / 2), label = my.dtw$index1[a],
       cex = 1.5)
  
  text(x = ((ref.dx[coord.ref[1]] + ref.dx[coord.ref[2]]) / 2), 
       y = ((ref.dz[coord.ref[1]] + ref.dz[coord.ref[2]]) / 2), label = my.dtw$index2[a],
       cex = 1.5)
}

my.dtw$costMatrix
my.dtw$normalizedDistance
plot(my.dtw)
plot(my.dtw,type="twoway")
dtwPlotDensity(my.dtw)

sqrt(sum((ml1[1,]-ml2[1,])^2))
