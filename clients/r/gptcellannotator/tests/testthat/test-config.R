test_that("gptca_config uses env defaults", {
  withr::with_envvar(
    c("GPTCA_CLI_PATH" = "/usr/local/bin/gca"),
    {
      cfg <- gptca_config()
      expect_s3_class(cfg, "GptcaConfig")
      expect_equal(cfg$cli_path, "/usr/local/bin/gca")
    }
  )
})

test_that("gptca_config_set stores config", {
  cfg <- gptca_config(cli_path = "/usr/bin/gca")
  prev <- gptca_config_set(cfg)
  expect_s3_class(prev, "GptcaConfig")
  expect_identical(gptca_config_get(), cfg)
})
