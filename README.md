# Central Dogma Explorer Shiny App

This Shiny app allows users to explore DNA transcription and translation with validation and interactive UI.

## Deployment to shinyapps.io

To deploy this app to shinyapps.io, follow these steps:

1. Install the `rsconnect` package if you haven't already:

```R
install.packages("rsconnect")
```

2. Load the package and set your account info (replace with your own token and secret):

```R
library(rsconnect)
rsconnect::setAccountInfo(name='your_account_name',
                          token='your_token',
                          secret='your_secret')
```

3. Deploy the app from the `translation` directory:

```R
rsconnect::deployApp('translation')
```

Make sure you have all dependencies installed and listed in the `DESCRIPTION` file.

## Dependencies

- shiny
- shinyjs
- rsconnect

## Notes

- The app is contained in `app.R` inside the `translation` directory.
- The `DESCRIPTION` file lists the package dependencies.
