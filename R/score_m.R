#' Score combination
#'
#' This function combines the survival distance score (scoreS) and the clinical independence score (scoreC) with weighted standardization.
#' @param scoreS Vector of survival distance scores.
#' @param scoreC Vector of clinical independence scores.
#' @param weight Weight to apply to score 2.
#' @keywords feature selection
#' @return Vector of scores with attribute names as names. Higher is better.
#' @export
#' @examples
#' scoreCombination()
scoreCombination = function(scoreS,scoreC,weight){
  #print(combo)
  # normalizing
  mu = mean(scoreS)
  std = sqrt(var(scoreS))
  scoreS_normal = (scoreS - mu)/std
  names(scoreS_normal) = names(scoreS)

  mu = mean(scoreC[scoreC != 0])
  std = sqrt(var(scoreC[scoreC != 0]))
  scoreC_normal = (scoreC - mu)/std
  names(scoreC_normal) = names(scoreC)

  #top 50 genes filtered by correlation
  sum_normal = scoreS_normal + scoreC_normal*weight
  names(sum_normal) = names(scoreS)

  return(sum_normal)
}
