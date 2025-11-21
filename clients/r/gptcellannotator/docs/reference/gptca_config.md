# Configure GPT Cell Annotator backend access

Create or manage configuration objects used to communicate with the GPT
Cell Annotator backend. The configuration stores base URLs (for optional
REST calls), authentication, CLI paths, and retry policies.

## Usage

``` r
gptca_config(
  base_url = Sys.getenv("GPTCA_BASE_URL", unset = NA_character_),
  api_key = Sys.getenv("GPTCA_API_KEY", unset = NA_character_),
  cli_path = Sys.getenv("GPTCA_CLI_PATH", unset = ""),
  timeout = 120,
  retry_max = 3,
  retry_backoff = 1,
  offline = TRUE,
  user_agent = utils::packageName()
)

gptca_config_set(config)

gptca_config_get(default = gptca_config())

gptca_config_reset()
```

## Arguments

- base_url:

  Optional base URL for the REST backend. Set to `NA` to skip HTTP
  calls.

- api_key:

  Optional API key passed as a Bearer token.

- cli_path:

  Path to the `gca` CLI (used by default); detected automatically when
  not supplied.

- timeout:

  Request timeout in seconds.

- retry_max:

  Maximum retry attempts for transient HTTP failures.

- retry_backoff:

  Initial delay in seconds applied to exponential backoff.

- offline:

  Logical flag forcing CLI mode (default `TRUE`).

- user_agent:

  Additional user-agent suffix appended to requests.

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
