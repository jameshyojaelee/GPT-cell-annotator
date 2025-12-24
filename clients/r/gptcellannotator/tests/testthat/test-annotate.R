fake_validated_response <- function() {
  list(
    summary = list(
      total_clusters = 1L,
      supported_clusters = 1L,
      flagged_clusters = 0L,
      unknown_clusters = list()
    ),
    clusters = list(
      list(
        cluster_id = "0",
        status = "supported",
        confidence = "High",
        warnings = list(),
        annotation = list(
          primary_label = "B cell",
          ontology_id = "CL:0000236",
          confidence = "High",
          rationale = "Markers support B cell identity",
          markers = list("MS4A1", "CD79A")
        ),
        validation = list(
          cluster_id = "0",
          primary_label = "B cell",
          is_supported = TRUE
        )
      )
    )
  )
}

test_that("gptca_prepare_clusters handles named lists", {
  clusters <- gptca_prepare_clusters(list(`0` = c("MS4A1", "CD79A")))
  expect_equal(length(clusters), 1L)
  expect_equal(clusters[[1]]$cluster_id, "0")
  expect_equal(clusters[[1]]$markers, c("MS4A1", "CD79A"))
})

test_that("gptca_prepare_clusters handles data frames", {
  df <- data.frame(
    cluster = c(0, 0, 1, 1),
    gene = c("MS4A1", "CD79A", "CD3E", "LCK"),
    avg_log2FC = c(2, 1.5, 1.8, 1.2)
  )
  clusters <- gptca_prepare_clusters(df, top_n = 1)
  expect_equal(length(clusters), 2L)
  expect_equal(clusters[[1]]$markers, "MS4A1")
})

test_that("gptca_annotate_markers uses CLI", {
  cfg <- gptca_config(cli_path = "/mock/cli", offline = TRUE)
  result <- with_mocked_bindings(
    gptca_annotate_markers(list(`0` = c("MS4A1", "CD79A")), config = cfg),
    gptca_cli_available = function(config) TRUE,
    gptca_cli_annotate = function(clusters, dataset_context, config, offline = TRUE, marker_limit = NULL) {
      fake_validated_response()
    }
  )
  expect_s3_class(result, "gptca_annotation")
  expect_true(result$validated)
  expect_equal(result$clusters$primary_label, "B cell")
  expect_equal(result$summary$total_clusters, 1L)
})
