#' Score s
#'
#' This function calculates the survival distance score based on attribute values and survival data.
#' @param attrib Matrix of atrribute values with attribute names as row names and sample names as column names.
#' @param surv data.frame of survival invormation with sample names as row names and column names 'overall.survival.indicator' and 'OVERALL.SURVIVAL'.
#' @keywords feature selection
#' @return Vector of scores with attribute names as names. Higher is better.
#' @export
#' @examples
#' get_score_s()
get_score_s <- function(attrib, surv){

  RowVar <- function(x, ...) {
    rowSums((x - rowMeans(x, ...))^2, ...)/(dim(x)[2] - 1)
  }
  sigma2 <- RowVar(attrib)
  # sorted timepoints of events
  timepoints <- sort(unique(surv[surv$overall.survival.indicator==1,]$OVERALL.SURVIVAL))
  # vector of weights for each gene
  sds <- rep(0, each=dim(attrib)[1])

  for (t in timepoints){
    m <- rowMeans(subset(attrib, select = rownames(surv[surv$OVERALL.SURVIVAL >= t,])))
    m[is.nan(m)] = 0 # this transforms all NaN into 0
    n <- dim(subset(attrib, select = rownames(surv[surv$OVERALL.SURVIVAL >= t,])))[2]


    mh <- rowMeans(subset(attrib, select = rownames(surv[surv$OVERALL.SURVIVAL <= t & surv$overall.survival.indicator==1,])))
    mh[is.nan(mh)] <- 0 # this transforms all NaN into 0
    nh <- dim(subset(attrib, select = rownames(surv[surv$OVERALL.SURVIVAL <= t & surv$overall.survival.indicator==1,])))[2]

    # function
    # look only at patients with event at time t
    patients <- rownames(surv[surv$OVERALL.SURVIVAL == t & surv$overall.survival.indicator == 1,])

    w1 <- rep(0, each=dim(attrib)[1])
    for (p in patients){
      w1 <- w1 + n*(attrib[,p] - m)^2 + nh*(attrib[,p] - mh)^2
    }
    w1 <- w1/length(patients)
    sds <- sds + w1

  }
  sds <- sds/sigma2
  return(sds)
}
