% Generated manually until roxygen can be run
\name{gptca_config}
\alias{gptca_config}
\alias{gptca_config_set}
\alias{gptca_config_get}
\alias{gptca_config_reset}
\title{Configure GPT Cell Annotator CLI access}
\description{
Create or manage configuration objects used to communicate with the GPT Cell Annotator CLI. The configuration stores the CLI path and offline flag.
}
\usage{
gptca_config(
  cli_path = Sys.getenv("GPTCA_CLI_PATH", unset = ""),
  offline = TRUE
)

gptca_config_set(config)

gptca_config_get(default = gptca_config())

gptca_config_reset()
}
\arguments{
\item{cli_path}{Path to the \code{gca} CLI (used by default); detected automatically when not supplied.}

\item{offline}{Logical flag forcing CLI mode (default \code{TRUE}).}

\item{config}{A \verb{GptcaConfig} object to activate.}

\item{default}{Default configuration returned when none is active.}
}
\value{
\code{gptca_config()} returns an object of class \verb{GptcaConfig}. \code{gptca_config_set()} invisibly returns the previous config. \code{gptca_config_get()} returns the currently active config. \code{gptca_config_reset()} resets to defaults and returns the new config.
}
\examples{
if (interactive()) {
  cfg <- gptca_config(cli_path = Sys.which("gca"))
  gptca_config_set(cfg)
}
}
