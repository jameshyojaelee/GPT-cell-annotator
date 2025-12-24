#' GPT Cell Annotator configuration
#'
#' @description
#' Creates or retrieves configuration used by `gptcellannotator` to talk to the
#' shell out to the CLI (default).
#'
#' @param cli_path Optional path to the CLI (`gca`). If `NULL`, the CLI is
#'   auto-detected using `Sys.which("gca")`.
#' @param offline Logical flag forcing the CLI path. Defaults to `TRUE`.
#'
#' @return A `GptcaConfig` object.
#' @export
gptca_config <- function(
  cli_path = Sys.getenv("GPTCA_CLI_PATH", unset = ""),
  offline = TRUE
) {
  if (!nzchar(cli_path)) {
    cli_path <- Sys.which("gca")
    if (!nzchar(cli_path)) {
      cli_path <- NA_character_
    }
  }

  cfg <- list(
    cli_path = if (!is.null(cli_path) && nzchar(cli_path)) cli_path else NULL,
    offline = isTRUE(offline)
  )
  class(cfg) <- "GptcaConfig"
  cfg
}

#' @export
print.GptcaConfig <- function(x, ...) {
  cli::cli_text("{.strong GPT Cell Annotator configuration}")
  cli::cli_ul(c(
    "cli_path" = ifelse(is.null(x$cli_path), "<not-found>", x$cli_path),
    "offline" = ifelse(isTRUE(x$offline), "yes", "no")
  ))
  invisible(x)
}

#' Activate a configuration
#'
#' @param config A `GptcaConfig` object.
#' @return The previously active configuration (invisibly).
#' @export
gptca_config_set <- function(config) {
  if (!inherits(config, "GptcaConfig")) {
    cli::cli_abort("{.arg config} must be a {.cls GptcaConfig}.")
  }
  prev <- gptca_config_get(default = NULL)
  assign("gptca_config", config, envir = .gptca_state)
  invisible(prev)
}

#' Retrieve the active configuration
#'
#' @param default Value returned if no config has been activated.
#' @return A `GptcaConfig` object.
#' @export
gptca_config_get <- function(default = gptca_config()) {
  cfg <- get0("gptca_config", envir = .gptca_state, ifnotfound = NULL)
  if (is.null(cfg)) {
    cfg <- default
    if (!is.null(cfg)) {
      gptca_config_set(cfg)
    }
  }
  cfg
}

#' Reset configuration to defaults
#'
#' @export
gptca_config_reset <- function() {
  assign("gptca_config", NULL, envir = .gptca_state)
  invisible(gptca_config_get())
}

.gptca_state <- new.env(parent = emptyenv())

`%||%` <- function(x, y) if (is.null(x)) y else x
