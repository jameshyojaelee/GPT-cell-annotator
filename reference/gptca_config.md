# Configure GPT Cell Annotator CLI access

Create or manage configuration objects used to communicate with the GPT
Cell Annotator CLI. The configuration stores the CLI path and offline
flag.

## Usage

``` r
gptca_config(
  cli_path = Sys.getenv("GPTCA_CLI_PATH", unset = ""),
  offline = TRUE
)

gptca_config_set(config)

gptca_config_get(default = gptca_config())

gptca_config_reset()
```

## Arguments

- cli_path:

  Path to the `gca` CLI (used by default); detected automatically when
  not supplied.

- offline:

  Logical flag forcing CLI mode (default `TRUE`).

- config:

  A `GptcaConfig` object to activate.

- default:

  Default configuration returned when none is active.

## Value

`gptca_config()` returns an object of class `GptcaConfig`.
`gptca_config_set()` invisibly returns the previous config.
`gptca_config_get()` returns the currently active config.
`gptca_config_reset()` resets to defaults and returns the new config.

## Examples

``` r
if (interactive()) {
  cfg <- gptca_config(cli_path = Sys.which("gca"))
  gptca_config_set(cfg)
}
```
