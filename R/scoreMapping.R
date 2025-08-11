#' Mapping scores
#'
#' @inheritParams plotSameLayout
#'
#' @return An S4Vectors list containing the cluster accuracy `ClusterAccuracy`,
#' the label precision (`LabelPrecision`), the stage precision (`StagePrecision`),
#' the global accuracy (`GlobalAccuracy`), the global precision (`GlobalPrecision`),
#' the mapping score (`MappingScore`)
#' @export
#'
#' @examples
scoreMapping <- function(net,
                         new_cor,
                         label_attr = 'SubClass') {

  if(is.character(label_attr) && length(label_attr) == 1 && label_attr %in% colnames(SingleCellExperiment::colData(net))) {
    label_attr <- SingleCellExperiment::colData(net)[,label_attr]
  }
  label_attr <- as.vector(label_attr)

  mean_cors <- t(get_column_group_average(t(new_cor), interaction(label_attr, net$Stages, sep = 'and')))

  accuracy <- apply(mean_cors, 2, max)

  precision <- apply(mean_cors, 1, max)

  max_cor <- stats::cor(t(apply(as.matrix(SingleCellExperiment::logcounts(net)[rownames(net)[SingleCellExperiment::rowData(net)$informative],]), 1, function(v) {
    (v-mean(v))/stats::sd(v)
  })), as.matrix(SingleCellExperiment::logcounts(net)[rownames(net)[SingleCellExperiment::rowData(net)$informative],]))

  max_mean_cors <- t(get_column_group_average(t(max_cor), interaction(label_attr, net$Stages, sep = 'and')))

  max_accuracy <- mean(apply(max_mean_cors, 2, max))

  accuracy <- accuracy/max_accuracy

  max_precision <- apply(max_mean_cors, 1, max)

  precision <- precision/max_precision

  by_subclass_precision <- tapply(precision, unlist(lapply(strsplit(names(precision), 'and', fixed = TRUE), function(i) {i[1]})), mean)
  by_stage_precision <- tapply(precision, unlist(lapply(strsplit(names(precision), 'and', fixed = TRUE), function(i) {i[2]})), mean)

  global_accuracy <- mean(accuracy)

  accuracy[which(accuracy > 1)] <- 1

  global_precision <- mean(precision)

  mapping_score <- mean(c(global_accuracy, global_precision))

  return(S4Vectors::List(ClusterAccuracy = accuracy,
                         LabelPrecision = sort(by_subclass_precision),
                         StagePrecision = sort(by_stage_precision),
                         GlobalAccuracy = global_accuracy,
                         GlobalPrecision = global_precision,
                         MappingScore = mapping_score))

}
