# Graphical user interface to application modules
# nocov start

ADSSModApp = function() {

  appDir = system.file("ADSSModApp", package = "MedianaDesigner")
  shiny::runApp(appDir, display.mode = "normal")

}

ADTreatSelApp = function() {

  appDir = system.file("ADTreatSelApp", package = "MedianaDesigner")
  shiny::runApp(appDir, display.mode = "normal")

}

ADPopSelApp = function() {

  appDir = system.file("ADPopSelApp", package = "MedianaDesigner")
  shiny::runApp(appDir, display.mode = "normal")

}

FutRuleApp = function() {

  appDir = system.file("FutRuleApp", package = "MedianaDesigner")
  shiny::runApp(appDir, display.mode = "normal")

}

EventPredApp = function() {

  appDir = system.file("EventPredApp", package = "MedianaDesigner")
  shiny::runApp(appDir, display.mode = "normal")

}

ADRandApp = function() {

  appDir = system.file("ADRandApp", package = "MedianaDesigner")
  shiny::runApp(appDir, display.mode = "normal")

}

MultAdjApp = function() {

  appDir = system.file("MultAdjApp", package = "MedianaDesigner")
  shiny::runApp(appDir, display.mode = "normal")

}

# nocov end
