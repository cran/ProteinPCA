#' @title Principal Component Analysis (PCA) Tool on Protein Expression Data
#' @name pca_analysis
#' @description This function performs PCA on protein expression data and produces a biplot using ggplot2.
#' Optionally, it supports grouping/clustering visualization with PCA loadings and confidence ellipses.
#'
#' @param data A numeric matrix or data frame of protein expression data. Rows are samples, columns are features.
#' @param scale Logical. Should the data be scaled? Default is TRUE.
#' @param center Logical. Should the data be centered? Default is TRUE.
#' @param plot Logical. Should a PCA biplot be generated? Default is TRUE.
#' @param groups Optional. A factor or character vector specifying group memberships of the samples.
#'
#' @return A list containing:
#' \item{pca}{The PCA object from \code{prcomp}.}
#' \item{explained_variance}{The percentage of variance explained by each principal component.}
#' \item{plot}{A ggplot2 PCA biplot (if \code{plot = TRUE}).}
#'
#' @examples
#' set.seed(123)
#' data_matrix <- matrix(rnorm(100 * 20), nrow = 100, ncol = 20)
#' rownames(data_matrix) <- paste0('Sample_', 1:100)
#' colnames(data_matrix) <- paste0('Protein_', 1:20)
#' groups <- sample(c("Group A", "Group B"), 100, replace = TRUE)
#' result <- pca_analysis(data_matrix, groups = groups)
#' print(result$explained_variance)
#' print(result$plot)
#'
#' @author Paul Angelo C. Manlapaz
#' @references
#' Jolliffe, I. (2001). Principal Component Analysis (2nd ed.). Springer. https://doi.org/10.1007/b98835
#' Gabriel, K. R. (1971). The biplot graphic display of matrices with application to principal component analysis. Biometrika, 58(3), 453–467. https://doi.org/10.1093/biomet/58.3.453
#' Zhang, Z., Chen, L., Sun, B., Ruan, Z., Pan, P., Zhang, W., Jiang, X., Zheng, S., Cheng, S., Xian, L., Wang, B., Yang, J., Zhang, B., Xu, P., Zhong, Z., Cheng, L., Ni, H., & Hong, Y. (2024). Identifying septic shock subgroups to tailor fluid strategies through multi-omics integration. Nature Communications, 15(1). https://doi.org/10.1038/s41467-024-53239-9
#' Anandan, A., Nagireddy, R., Sabarinathan, S., Bhatta, B. B., Mahender, A., Vinothkumar, M., Parameswaran, C., Panneerselvam, P., Subudhi, H., Meher, J., Bose, L. K., & Ali, J. (2022). Multi-trait association study identifies loci associated with tolerance of low phosphorus in Oryza sativa and its wild relatives. Scientific Reports, 12(1). https://doi.org/10.1038/s41598-022-07781-5
#'
#' @export
#' @importFrom stats prcomp
#' @importFrom ggplot2 ggplot aes geom_point geom_text ggtitle xlab ylab stat_ellipse theme_minimal theme element_text margin scale_color_brewer geom_segment geom_label
#' @importFrom grid arrow unit
#' @importFrom gridExtra grid.arrange

utils::globalVariables(c("PC", "PC1", "PC2", "Samples", "Groups", "xend", "yend", "Variance", "Variables"))

pca_analysis <- function(data, scale = TRUE, center = TRUE, plot = TRUE, groups = NULL) {

  # Validate input
  if (!is.data.frame(data) && !is.matrix(data)) {
    stop("Input data must be a data frame or matrix.")
  }

  # Run PCA
  pca_res <- stats::prcomp(data, scale. = scale, center = center)
  explained_var <- (pca_res$sdev^2) / sum(pca_res$sdev^2) * 100
  result <- list(pca = pca_res, explained_variance = explained_var)

  if (plot) {
    # Check for required packages
    if (!requireNamespace("ggplot2", quietly = TRUE) || !requireNamespace("gridExtra", quietly = TRUE)) {
      stop("Both ggplot2 and gridExtra are required for plotting. Please install them with install.packages().")
    }

    # Prepare PCA scores
    df_scores <- as.data.frame(pca_res$x[, 1:2])
    df_scores$Samples <- rownames(df_scores)

    # Handle grouping
    if (!is.null(groups)) {
      if (length(groups) != nrow(df_scores)) stop("Length of 'groups' must match number of samples.")
      df_scores$Groups <- as.factor(groups)
    }

    # PCA loadings
    loadings <- as.data.frame(pca_res$rotation[, 1:2])
    loadings$Variables <- rownames(loadings)

    # Scale loadings for arrow plotting
    scale_factor <- max(abs(df_scores$PC1), abs(df_scores$PC2)) * 0.8
    loadings$xend <- loadings$PC1 * scale_factor
    loadings$yend <- loadings$PC2 * scale_factor

    # --- Biplot ---
    p1 <- ggplot2::ggplot(df_scores, ggplot2::aes(x = PC1, y = PC2)) +
      ggplot2::xlab(paste0("PC1 (", round(explained_var[1], 2), "%)")) +
      ggplot2::ylab(paste0("PC2 (", round(explained_var[2], 2), "%)")) +
      ggplot2::ggtitle("PCA Biplot of Protein Expression Data") +
      ggplot2::theme_minimal(base_size = 14) +
      ggplot2::theme(
        plot.title = ggplot2::element_text(face = "bold", size = 20, hjust = 0.5, margin = ggplot2::margin(b = 15)),
        axis.title = ggplot2::element_text(face = "bold", size = 16, margin = ggplot2::margin(t = 10, r = 10))
      )

    # Plot sample points and ellipses if groups exist
    if (!is.null(groups)) {
      p1 <- p1 +
        ggplot2::geom_point(ggplot2::aes(color = Groups), size = 3) +
        ggplot2::stat_ellipse(ggplot2::aes(color = Groups), level = 0.95, linetype = "solid", linewidth = 1, alpha = 0.2) +
        ggplot2::scale_color_brewer(palette = "Set1")
    } else {
      p1 <- p1 + ggplot2::geom_point(color = "#766CDB", size = 3)
    }

    # Sample labels
    p1 <- p1 + ggplot2::geom_text(ggplot2::aes(label = Samples), vjust = -0.5, size = 3)

    # Add loadings arrows and labels
    p1 <- p1 +
      ggplot2::geom_segment(data = loadings,
                            ggplot2::aes(x = 0, y = 0, xend = xend, yend = yend),
                            arrow = grid::arrow(length = grid::unit(0.3, "cm")),
                            color = "gray40", linewidth = 0.8) +
      ggplot2::geom_label(data = loadings,
                          ggplot2::aes(x = xend, y = yend, label = Variables),
                          color = "gray10", fill = "white", fontface = "bold", size = 3)

    # --- Scree Plot (Top 10 PCs) ---
    df_var <- data.frame(PC = paste0("PC", seq_along(explained_var)), Variance = explained_var)
    p2 <- ggplot2::ggplot(df_var[1:10, ], ggplot2::aes(x = PC, y = Variance)) +
      ggplot2::geom_bar(stat = "identity", fill = "#69b3a2") +
      ggplot2::theme_minimal() +
      ggplot2::ggtitle("Scree Plot (Top 10 PCs)") +
      ggplot2::xlab("") +
      ggplot2::ylab("Variance Explained (%)") +
      ggplot2::theme(
        axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
        plot.title = ggplot2::element_text(face = "bold", hjust = 0.5)
      )

    # Combine plots
    result$plot <- gridExtra::grid.arrange(p1, p2, ncol = 2)
  }

  return(result)
}
