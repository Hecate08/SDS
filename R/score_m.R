#' Score combination
#'
#' This function combines the survival distance score (score1) and the clinical independence score (score2) with weighted standardization.
#' @param score1 Vector of survival distance scores.
#' @param score2 Vector of clinical independence scores.
#' @param weight Weight to apply to score 2.
#' @keywords feature selection
#' @return Vector of scores with attribute names as names. Higher is better.
#' @export
#' @examples
#' scoreCombination()
scoreCombination = function(score1,score2,weight){
  #print(combo)
  # normalizing
  mu = mean(score1)
  std = sqrt(var(score1))
  score1_normal = (score1 - mu)/std
  names(score1_normal) = names(score1)
  #plot(density(score1_normal))

  mu = mean(score2[score2 != 0])
  std = sqrt(var(score2[score2 != 0]))
  score2_normal = (score2 - mu)/std
  names(score2_normal) = names(score2)

  #top 50 genes filtered by correlation
  sum_normal = score1_normal + score2_normal*weight
  names(sum_normal) = names(score1)

  return(sum_normal)
}
