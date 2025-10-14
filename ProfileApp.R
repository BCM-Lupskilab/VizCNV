library(profvis)
app_dir <- "VizCNV"   # path to your app
profvis({
  shiny::runApp(appDir = app_dir, host = "127.0.0.1", port = 5050, launch.browser = T)
})
