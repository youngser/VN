suppressMessages(require(VN))
suppressMessages(require(igraph))
suppressMessages(require(ggplot2))
suppressMessages(require(reshape2))
suppressMessages(require(dplyr))

suppressMessages(require(parallel))
suppressMessages(require(doMC))

registerDoMC(cores=detectCores()-1)
getDoParWorkers()


doMC <- function(corr=0.9,sim=TRUE,HS.full=FALSE,HA=FALSE)
{
    s <- 4 # number of seeds for SGM, assuming that they are the first s vertices!
    x <- 1 # voi in G1, to be found in G2
    h <- ell <-  1 # max walk from voi -- all vertices in nbd should be within path of length h from at least one seed
    gamma <- 0.1 # gamma <- max tolerance for alpha, how far away from the barycenter user is willing to go for
         # the initialization of SGM on any given iteration

    R <- 100 # number of times to run sgm
    numMC <- 1000 # number of ER
    nc <- 1#16

    i <- 1
    mc <- 0
    require(VN)

    s <- ifelse(sim,4,12)

    m <- 5 # junk on G1
    n <- 20 # |V| = 1 + s + n + m = 25
    (nV1 <- 1+s+n+m)
    mp <- 0 # junk on G2, |V'| = 1 + s + n + m', remove m from G2
    (nV2 <- 1+s+n+mp)
    d <- 5
    p <- 0.5

    MC <- foreach(mc=1:numMC) %dopar% {
### generate a pair of correlated graphs
#set.seed(246) # n=10
#set.seed(146+mc) # n=20, html mc=0
        set.seed(32456*corr+mc) # n=20
#    set.seed(1346*corr+mc) # n=20
#    gg <- sample_correlated_gnp_pair(nV1, corr, p, directed = FALSE)
        if (sim) {
            lpvs <- sample_sphere_surface(dim=d, n=nV1)/1.5
            gg <- rdpg.sample.correlated(t(lpvs),corr)
#    ggg <- rg.sample.SBM.correlated(nV1, B=matrix(c(0.7,0.1,0.1,0.7), nrow=2), rho=c(0.3,0.7), sigma=corr,conditional=TRUE)
#    gg <- ggg[[1]]

            g1 <- gg[[1]]; #summary(g1)
            g2 <- gg[[2]]; #summary(g2)
#    g1 <- sample_dot_product(lpvs) # rdpg
#    g2 <- sample_correlated_gnp(g1, corr=corr)
#    cor(as.vector(g1[]), as.vector(g2[]))
        ## make g2 smaller than g1
            g2 <- delete_vertices(g2,v=(nV2+1):nV1); #summary(g2)
            W <- intersect(igraph::V(g1),igraph::V(g2))
        } else {
            data("HSgraphs")
#        print(load("~/Dropbox/SGM/nbdmaking/HSgraphfull.Rbin"))
            if (HS.full) {
                perm.fb <- c(coremap[,1],setdiff(1:vcount(HSfbgraphfull),coremap[,1]))
                perm.fr <- c(coremap[,2],setdiff(1:vcount(HSfriendsgraphfull),coremap[,2]))
                g1 <- permute.vertices(HSfbgraphfull,perm.fb) # (156,1437)
                g2 <- permute.vertices(HSfriendsgraphfull,perm.fr) # (134,668)
                W <- 1:nrow(coremap) # sees should be selected from the core (shared & corresponding) vertices
            } else {
                g1 <- HSfbgraphcore
                g2 <- HSfrgraphcore
                W <- intersect(igraph::V(g1),igraph::V(g2))
            }
        }

    # pick VOI from the smaller graph (shared, W)
        (x <- sample(W,1))
        W <- setdiff(W,x) # exclude x from W

    # choose s seeds from W
        (S <- getSeeds(W,x,s))

## ----sgm,echo=TRUE,eval=TRUE---------------------------------------------
        NBDS <- vnsgm.ordered(x,S,g1,g2,h,ell,R,gamma,verb=FALSE)
        (case <- NBDS$case)
        seed <- NBDS$Sx #NBDS$labelsGx[1:s]
        s <- length(seed)

        if (case=="possible") {
            cat("mc =", mc, " is ", case, "!\n")

            vstar <- NBDS$labelsGx[s+1]
            vstar.ind <- which(NBDS$labelsGx==x)
#    candidate <- vstar.ind:length(NBDS$labelsGxp) # NBDS$labelsGxp[-c(1:vstar)]
            candidate <- sort(unique(c(vstar.ind,match(NBDS$Cxp,NBDS$labelsGxp)))) # x.ind:length(NBDS$labelsGxp) # NBDS$labelsGxp[-c(1:vstar)]
            ncandidate <- length(candidate)

#            Nx <- sort(unlist(ego(g1,h,nodes=x,mindist=0)))
#            Nxp <- sort(unlist(ego(g2,h,nodes=x,mindist=0)))
#            Nxp2 <- intersect(Nx,Nxp)
#            Exp <- sum(g2[][x,seed])

            prob <- NBDS$P[vstar.ind,candidate]
            names(prob) <- NBDS$labelsGxp[candidate]
    #prob
            vstar <- prob[which.max(prob)]
#vstar
#plot(as.integer(names(prob)), prob, type="h",col=2, lwd=2)
            vstar <- as.integer(names(vstar))

            rank.prob <- rank(-prob,ties.method = "average")
            rank.vstar <- rank.prob[1]

            ## if HA==TRUE, shuffle the prob and calc nrank
            if (HA) {
                rank.vstar2 <- rep(0,numMC)
                for (i in 1:numMC) {
                    prob2 <- sample(prob)
                    rank.prob2 <- rank(-prob2,ties.method = "average")
                    rank.vstar2[i] <- rank.prob2[1]
                }
                nrank <- (rank.vstar-1) / (ncandidate-1)
                nrank2 <- (rank.vstar2-1) / (ncandidate-1)
                (pval <- sum(nrank > nrank2) / numMC)
            } else {
                pval <- 0
            }
        } else { # impossible1
            seed <- Nxp2 <- Exp <- vstar <- rank.vstar <- pval <- NA
            if (case=="impossible2") {
                ncandidate <- length(NBDS$Cxp)
            } else {
                ncandidate <- 0
            }
        }
        cat("done mc=", mc, "\n")

        c(case, length(seed), ncandidate, rank.vstar, pval, vstar)#, length(Nxp2), Exp)
    }

    return(MC)
}


doSim <- function()
{
    system.time(MC <- doMC(0.1, HA=TRUE)); MC1 <- MC
#    save(MC, file="MC-RDPG-rho0.1-seed1-N1-new2.Rbin")
    system.time(MC <- doMC(0.9, HA=TRUE)); MC9 <- MC
#    save(MC, file="MC-RDPG-rho0.9-seed1-N1-new2.Rbin")

#    load("MC-RDPG-rho0.1-seed1-N1-new2.Rbin"); MC1 <- MC
#    load("MC-RDPG-rho0.9-seed1-N1-new2.Rbin"); MC9 <- MC
    case1 <- (sapply(MC1,"[",1))
    case9 <- (sapply(MC9,"[",1))
    seed1 <- paste0("s",sapply(MC1,"[",2)); seed1[case1!="possible"] <- "s0"
    seed9 <- paste0("s",sapply(MC9,"[",2)); seed9[case9!="possible"] <- "s0"
    ncandidate1 <- as.numeric(sapply(MC1,"[",3))
    ncandidate9 <- as.numeric(sapply(MC9,"[",3))
    rank.vstar1 <- as.numeric(sapply(MC1,"[",4))
    rank.vstar9 <- as.numeric(sapply(MC9,"[",4))
    df1 <- data.frame(corr="0.1",case=case1,seed=seed1,ncand=ncandidate1,rank=rank.vstar1)
    df1 <- df1 %>% mutate(norm.rank=(rank-1)/(ncand-1))
    df9 <- data.frame(corr="0.9",case=case9,seed=seed9,ncand=ncandidate9,rank=rank.vstar9)
    df9 <- df9 %>% mutate(norm.rank=(rank-1)/(ncand-1))
    df19 <- rbind(df1,df9)
    df191 <- subset(df19,case=="possible")
#    kable(df1 %>% group_by(case,seed) %>% summarize(count=n(),mean.ncand=mean(ncand), mean.nrank=mean(norm.rank)),caption="A summary of the simulation for `corr=0.1`")
#    kable(df9 %>% group_by(case,seed) %>% summarize(count=n(),mean.ncand=mean(ncand), mean.nrank=mean(norm.rank,na.rm=TRUE)),caption="A summary of the simulation for `corr=0.9`")

    ggplot(df191,aes(x=ncand,fill=corr)) + scale_fill_brewer(palette="Set1") +
        geom_histogram(binwidth=1,position="identity",alpha=0.5)
    ggplot(df191,aes(x=ncand,fill=corr)) + scale_fill_brewer(palette="Set1") +
        facet_wrap(~seed) +
        geom_histogram(binwidth=1,position="identity",alpha=0.5)

    df191$ncand <- paste0("c",sprintf("%02d",df191$ncand))
    df191 <- subset(df191,case=="possible")

    ggplot(df191,aes(x=norm.rank,fill=corr)) + scale_fill_brewer(palette="Set1") +
        facet_wrap(~seed) +
        geom_histogram(binwidth=.05,position="identity",alpha=0.5)#, aes(fill=..count..))

    ggplot(subset(df191,seed %in% c("s1","s2","s3")),aes(x=norm.rank, fill=corr)) +
        scale_fill_brewer(palette="Set1") +
        facet_grid(ncand~seed) +
        geom_histogram(binwidth=.05,position="identity",alpha=0.5)#, aes(fill=..count..))
}

