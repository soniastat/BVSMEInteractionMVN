
# This updates the namespace
devtools::document() # As we use roxygen2


devtools::check() # Fix any errors or warnings if possible


devtools::build() # Build the package

# Then from the Git, push to commit and change in the github


# install the package
devtools::install_github("soniastat/BVSMEInteractionMVN")
