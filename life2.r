#Mrandom <- function(M,N) { matrix (as.integer(2*rnorm(N*M))%%2,M,N) }
Mrandom <- function(M,N,p) { matrix (as.integer(rexp(N*M)<p),M,N) }

nbs <- function(M,i,j) 
{ 
 istart <- max(1,i-1)
 jstart <- max(1,j-1)
 
 iend <- min(i+1,nrow(M))
 jend <- min(j+1,ncol(M))
 H<-M[seq(istart,iend,1),seq(jstart,jend,1)];
 H
}

update <- function(j,i,G)
{
 z <- G[i,j]			# z is cell status, 0 for dead and 1 for live
 y <- sum(nbs(G,i,j))- z	# y is the number of live neighbours
 
 cell_status_next <- z*(as.integer(y<2 | y>3)*0 + as.integer(y==2 | y==3)*1) + (1-z)*(as.integer(y==3)*1 + as.integer(y!=3)*0)
 cell_status_next
}

update_j <- function(i,N,G) { sapply(1:N,update,i,G) }

update_ij <- function(G,M,N)
{
 R <-matrix(0,M,N)
 R <- sapply(1:M,update_j,N,G)
 t(R)
}

cell.display <- function(G,S)
{
 #Sys.sleep(0.5)
 png(S,width=1024,height=1024,bg='black');		# Set output as png
 par(new=F)
 image(1:ncol(G), 1:nrow(G), t(G), col=c('black','green'),axes=FALSE,xlab='',ylab='')
 dev.off()
}

life <- function (n_gen,M,N,p,wd='/tmp')
{
 set.seed(1337)
 this_gen <- Mrandom(M,N,p) 
 next_gen <- matrix(0,M,N)
 for (k in 1:n_gen)
 	{
	 print(k)
	 cell.display(this_gen,sprintf("%s/L%d.png",wd,k))
	 next_gen <- update_ij(this_gen,M,N)
	 this_gen <- next_gen
	}
}
