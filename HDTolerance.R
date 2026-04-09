## Univariate case
## Calculate a parameter needed to get Hellinger Distance for certian distribution
## For Normal z score returned
## For Poisson lambda is returned (assuming first lambda is 1)
## For beta, the value returned (x) means that the first distribution should use
## parameters 16/x and x and the second distribution x and 16/x
GetParametersForDistance<-function(d, distType)
{
    if (d<0 || d>1)
    {
        stop("Hellinger distance must be between 0 and 1")
    }
    switch(as.character(distType), 
           Normal=GetZForHellingerdistance(d),
           Poisson=GetLambdaForHellingerdistance(d),
           Beta=GetBetaParamsForHellingerdistance(d)
           )
}

# calculates z score between means assuming standard deviation is the same.
# generater vectors will have means M and M+z*SD, where M and SD are arbitarbly choosen
GetZForHellingerdistance<-function(d)
{
    if (d<0 || d>1)
    {
        stop("Hellinger distance must be between 0 and 1")
    }
    if (d==1)
    {
        return(10)
    }
    sqrt(-8*log(1-d*d))
}

# Calculates Lambda for Y vector, assuming lambda for X is 1
GetLambdaForHellingerdistance<-function(d)
{
    if (d<0 || d>1)
    {
        stop("Hellinger distance must be between 0 and 1")
    }
    if (d==1)
    {
        return(25)
    }
    (sqrt(-2*log(1-d*d))+1)^2
}

# returns alpha. Parameters for X are 16/a and a, for Y parameters are a and 16/a
GetBetaParamsForHellingerdistance<-function(d)
{
    if (d<0 || d>1)
    {
        stop("Hellinger distance must be between 0 and 1")
    }
    if (d==1)
    {
        return(1)
    }
    if (d==0)
    {
        return(4)
    }
    f<-function(a) TheoreticalHellingerBeta(16/a,a,a,16/a)-d
    alpha<-uniroot(f, c(1,4))
    alpha$root
}

#returns alpha. Second vector should use 1
GetExpParamsForGellingerDistance<-function(d)
{
    if (d<0 || d>1)
    {
        stop("Hellinger distance must be between 0 and 1")
    }
    if (d==1)
    {
        return(1e6)
    }
    b=2-(2/(1-d^2))^2
    -b/2+sqrt(b^2/4-1)
}

# Finds uniform distribution that yield desired HD
GetUnifParamsForHellingerDistance<-function(d, minX=0, maxX=1)
{
    if (d<0 || d>1)
    {
        stop("Hellinger distance must be between 0 and 1")
    }
    if (maxX<=minX)
    {
        stop("Right side of the first interval must be greater than the left side")
    }
    if (d==0)
    {
        return(c(0,1,0,1))
    }
    if (d==1)
    {
        return(c(1.1,2.1))
    }
    len1 = maxX-minX
    minY=d*len1*1.01+minX
    IntersectToGeomAverage = 1-d*d;
    intersect=max(0, maxX-minY)
    maxY = minY + (intersect/IntersectToGeomAverage)^2/(len1)
    list(FirstLeft=minX, FirstRight=maxX, SecondLeft=minY, SecondRight=maxY)
}

# Calculate theoretical HD for known univariate normal distributions
UnivariateNormalHellingerDistance<-function(mu1, mu2, sigma1, sigma2)
{
    totalVar = sigma1^2+sigma2^2
    sqrt(1-exp(-(mu1-mu2)^2/4/(totalVar))*sqrt(2*sigma1*sigma2/totalVar))
}

# Calculate theoretical HD for known univariate Poisson distributions
TheoreticalHellingerPoison<-function(L1,L2)
{
    sqrt(1-exp(-(L1+L2-2*sqrt(L1*L2))/2))
}

# Calculate theoretical HD for known univariate beta distributions
TheoreticalHellingerBeta<-function(Alpha1, Beta1, Alpha2, Beta2)
{
    sqrt(1-beta((Alpha1+Alpha2)/2, (Beta1+Beta2)/2)/sqrt(beta(Alpha1,Beta1)*beta(Alpha2,Beta2)))
}

# Calculate theoretical HD for known univariate exponential distributions
TheoreticalHellingerExponential<-function(L1,L2)
{
    sqrt(1-2*sqrt(L1*L2)/(L1+L2))
}

# Calculate theoretical HD for known uniform distributions
TheoreticalHellingerUniForm<-function(minX, maxX, minY, maxY)
{
    if (minX>=maxX || minY>=maxY)
    {
        stop("Invalid minimum and maximum values")
    }
    len1<-maxX-minX
    len2<-maxY-minY
    intersect = min(maxX, maxY) - max(minX, minY)
    if (intersect <= 0)
        return(1)
    sqrt(1-intersect/sqrt(len1*len2))
}

# given two densities estimated at different points, 
##approximate them at the other density points to get an alighed values
# Later used to calculate hellinger integral

AlignDensityFunctions<-function(x1, y1, x2, y2)
{
    One.Approximation<-approx(x1, y1, xout = x2, yleft = 0, yright = 0, rule = 2)
    Two.Approximation<-approx(x2, y2, xout = x1, yleft = 0, yright = 0, rule = 2)
    One.Order<-order(c(x1, One.Approximation$x))
    X1<-c(x1, One.Approximation$x)[One.Order]
    Y1<-c(y1, One.Approximation$y)[One.Order]
    Y2<-c(Two.Approximation$y, y2)[One.Order]
    list(x1=X1,y1=Y1,y2=Y2)
}

# given two alligned densities calculate hellinger integral
CalculateHellinger<-function(X, Y1, Y2)
{
    require(sfsmisc)
    n1 = length(X)
    n2 = length(Y1)
    if (n1 != n2)
        stop("Length of data point differs from length of PDF values")
    n3 = length(Y2)
    if (n2 != n3)
        stop("Densities have different number of points")
    combinedY<-sqrt(Y1*Y2)
    require(sfsmisc, quietly = TRUE)
    HellingerIntegral<-integrate.xy(X, combinedY, use.spline = FALSE)
    
    sqrt((1-min(HellingerIntegral,1)))
}

#Estimate hellinger distance between two univariate samples
UnivariateHellingerDistance<-function(X, Y, kernel="normal", type="KDE")
{
    if(type=="KDE")
    {
        # 1. Calculate densities
        require(KernSmooth)
        X.density<-bkde(X, kernel = kernel, gridsize = min(ceiling(length(X)/5), 400)+1)
        X.density$y[X.density$y<0]<-0
        Y.density<-bkde(Y, kernel = kernel, gridsize = min(ceiling(length(Y)/5), 400)+1)
        Y.density$y[Y.density$y<0]<-0
        # 2. Align densities
        AlignedDensities<-AlignDensityFunctions(X.density$x, X.density$y, Y.density$x, Y.density$y)
        # 3. Calculate hellinger integral and distance
        CalculateHellinger(AlignedDensities$x1, AlignedDensities$y1, AlignedDensities$y2)
    }
    else if (type=="hist")
    {
        H<-hist(c(X,Y), plot = FALSE, breaks=100)
        Points<-H$breaks
        Freq1<-hist(X, breaks = Points, plot = FALSE)$counts
        Freq2<-hist(Y, breaks = Points, plot = FALSE)$counts
        Freq1<-Freq1/sum(Freq1)
        Freq2<-Freq2/sum(Freq2)
        sqrt(1-sum(sqrt(Freq1*Freq2)))
    }
}

# Generate univariate random vectpr
generateRandomdistribtuion<-function(size, distributionType)
{
    switch(distributionType, 
           Normal = rnorm(size, mean = 5, sd=10),
           Poisson = rpois(size, lambda = 10),
           Beta = rbeta(size, shape1 = 2, shape2 = 5),
           Exponential = rexp(size, r=10),
           LogNormal =rlnorm(n=size, meanlog=5, sdlog=0.5),
           ChiSquared = rchisq(n=size, df=5))    
}

# Generate two column of data given size, distribution type 
# and a parameter that produces desired HD
generateRandomVectors<-function(size, distributionType, param)
{
    switch(as.character(distributionType), 
           Normal = data.frame(X=rnorm(size, mean = 5, sd=2),Y=rnorm(size, mean = 5+2*param, sd=2)),
           Poisson = data.frame(X=rpois(size, lambda = 1), Y=rpois(size, lambda = param)),
           Beta = data.frame(X=rbeta(size, shape1 = 16/param, shape2 = param), Y=rbeta(size, shape1 = param, shape2 = 16/param)),
           Exponential = data.frame(X=rexp(size, r=param), Y=rexp(size, r=1)))    
    
}

# Generate pairs of vectors that should have zero HD (coming from the same distribution)
# to analyse sensitive of algorithm to the sample size
GenerateDataForZeroDistance<-function(minSize=1000, maxSize=20000, step=500, SimNum=100, 
                       distributionTypes = c("Normal", "Poisson", "Beta", "Exponential", "LogNormal", "ChiSquared"))
{
    set.seed(1)
    
    
    frameSize=((maxSize-minSize)/step+1)*SimNum*length(distributionTypes)
    df<-expand.grid(SimSize=seq(minSize,maxSize,step),Run=1:SimNum,DistType=distributionTypes)
    params<-sapply(distributionTypes, function(x) GetParametersForDistance(0, x))
    for (distT in distributionTypes)
    {
        df$param[df$DistType==distT]=params[distT]
    }
    df$EmpiricalDistance<-mapply(GenerateAndCalculate, df$SimSize, df$DistType, df$param)
    df
}

# Generate pairs of vectors that should have certain HD 
# to test how close the algorithm to a "true" value
HellingerTheoreticalVersusEmpirical<-function(minSize=1000, maxSize=20000, step=500, SimNum=100, HDistances=seq(0,1,0.01),
                                              distributionTypes = c("Normal", "Poisson", "Beta"))

{
    set.seed(1)
    require(KernSmooth)
    require(dplyr)
    frameSize<-length(HDistances)*length(distributionTypes)*((maxSize-minSize)/step+1)*SimNum
    df<-expand.grid(SimSize=seq(minSize,maxSize,step),Run=1:SimNum,DistType=distributionTypes, TheoreticalDistance=HDistances)
    params<-expand.grid(DistType=distributionTypes, TheoreticalDistance=HDistances)
    params$param<-mapply(GetParametersForDistance, params$TheoreticalDistance, params$DistType)
    df<-left_join(df, params, by=c("TheoreticalDistance", "DistType"))
    df$EmpiricalDistance<-mapply(GenerateAndCalculate, df$SimSize, df$DistType, df$param, df$DistType=="Poisson")
    df
}

## Generate vector pair and estimate their HD
GenerateAndCalculate<-function(SimSize, DistType, param, useHist=FALSE)
{
    generatedVectors<-generateRandomVectors(SimSize, DistType, param)
    UnivariateHellingerDistance(generatedVectors$X, generatedVectors$Y)
    if (useHist)
    {
        UnivariateHellingerDistance(generatedVectors$X, generatedVectors$Y, type = "hist")
    }
    else
    {
        UnivariateHellingerDistance(generatedVectors$X, generatedVectors$Y)
    }
}

## Generate two datasets with the same number of columns from multivariate normal
## and estimate their HD using PCA decomposition. Used for calibration
GenerateAndCalculateMatrices<-function(SimSize, DimNo)
{
    covmat<-matrix(runif(DimNo^2),ncol=DimNo) #random non-negative matrix
    covmat<-crossprod(covmat) #make it symmetric and positive semi-definite
    e<-eigen(covmat)
    if(any(e$val<1E-6))
    {
        e$values[e$values<1E-6]<-1E-6 # ensure that matrix is positive definite
        CovMat<-e$vec %*% diag(e$val) %*% t(e$vec)
    }
    X<-GenerateMultiNormal(SimSize, DimNo, rep(0, DimNo), e)
    Y<-GenerateMultiNormal(SimSize, DimNo, rep(0, DimNo), e)
    HellingerPCADistance(X,Y)
}

## Generate data to calibrate HD
## data generated for different sample sizes
## HD is calculated for each PCA component
## Generated table is analysed to find the dependence between sample size and tolerance level
## (e.g. confidence interval for estimated HD where "True" HD must be zero )
GenerateZeroHDMultiNormal<-function(DimNo=10, SimSizes=c(seq(100,450,50),seq(500,900,100), seq(1000,1800,200),seq(2000,20000,500)), simNo=10)
{
    set.seed(1)
    df<-expand.grid(SimSize=SimSizes,Run=1:simNo)
    df<-cbind(df,t(mapply(GenerateAndCalculateMatrices, SimSize=df$SimSize, DimNo=DimNo)))
    df
    
}

## given a dataset calculate a matrix of HD between each pair of columns
CreateMatrixOfHD<-function(df)
{
    perm<-t(combn(1:ncol(df),2))
    v<-mapply(function(x,y) UnivariateHellingerDistance(df[,x], df[,y]), perm[,1], perm[,2])
    m<-diag(ncol(df))
    m[lower.tri(m)]<-v
    m[upper.tri(m)]<-t(m)[upper.tri(m)]
    m
}

## Fit regression to the list of zero HD results
## the regressionn is log(Y)~log(X) with heteroscedasticity
## weights are proportional to log(X)
AnalysisOfZeroHD<-function(SimSizes, MeasuredDistance, desiredQuantiles = c(P95=0.95,P99=0.99))
{
    logSimSize<-log(SimSizes)
    logMeaserDistance<-log(MeasuredDistance)
    X<-rep(logSimSize, ncol(logMeaserDistance)) #stack estimated distances against sample sizes
    Y<-stack(logMeaserDistance)[,1]
    
    lmfit<-lm(Y~X) #It is a good fit, but heteroscedasticity must be addressed
    # Plot of residuals vs predictors has a clear megaphone shape
    # Fitted values of regressing absolute values of residuals against the predictor is the estimate of sd
    lmWeights<-lm(abs(resid(lmfit))~X) #weight are resiprical of residual variances
    lmfitW<-lm(Y~X, weights = fitted.values(lmWeights)^-2)
    # get quantile values for max ZScore and chi Square
    predictors<-cbind(rep(1, length(logSimSize)), logSimSize)
    logMeans<-(predictors %*% coef(lmfitW))[,1]
    weights<-(predictors %*% coef(lmWeights))[,1]
    logerror<-sqrt(diag(predictors %*% vcov(lmfitW) %*% t(predictors))+summary(lmfitW)$sigma^2)*weights
    maxZScore<-apply((logMeaserDistance-logMeans)/logerror, 1, max)
    ChiSquares<-apply((logMeaserDistance-logMeans)/logerror, 1, function(x) sum(x[x>0]^2))
    ZScore=quantile(maxZScore, desiredQuantiles)
    if (length(names(desiredQuantiles))==length(ZScore))
    {
        names(ZScore)<-names(desiredQuantiles)
    }
    list(coef = coef(lmfitW), coefvar=vcov(lmfitW), se=summary(lmfitW)$sigma, weightCoef=coef(lmWeights), ZScore=ZScore, ChiSquare=quantile(ChiSquares, desiredQuantiles))
}

getLogMeanOfHDtolerance<-function(SimSizes, params)
{
    params$coef[1]+params$coef[2]*log(SimSizes)
}

getLogSDOfHDTolerance<-function(SimSizes, params)
{
    predictors<-cbind(rep(1, length(SimSizes)), log(SimSizes))
    weights<-(predictors %*% params$weightCoef)[,1]
    logerror<-sqrt(diag(predictors %*% params$coefvar %*% t(predictors))+params$se^2)*weights
    logerror
}

getLimitsOfHD<-function(SimSizes, params)
{
    logMeans<-getLogMeanOfHDtolerance(SimSizes, params) 
    logError<-getLogSDOfHDTolerance(SimSizes, params)
    retval<-logError %*% t(params$ZScore)
    retval<-sweep(retval, 1, logMeans, FUN = "+")
    exp(retval)
}

getMaxZScoreAndChisq<-function(vecDist, logMean, logError)
{
    ZScores<-(log(vecDist)-logMean)/logError
    list(max(ZScores), sum((ZScores[ZScores>0])^2))
}

createPlot<-function(df, simSize=-1)
{
    require(ggplot2)
    if (simSize<0)
    {
        g<-ggplot(data = df, aes(TheoreticalDistance, EmpiricalDistance-TheoreticalDistance))        
    }
    else
    {
        g<-ggplot(data = df[df$SimSize==simSize,], aes(TheoreticalDistance, EmpiricalDistance-TheoreticalDistance)) 
    }
    g<-g+geom_point()+facet_grid(DistType ~ .)
    g
}
savePlot<-function(file, object)
{
    png(file)
    object
    dev.off()
}
