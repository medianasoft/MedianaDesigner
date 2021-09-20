require(devEMF)
require(officer)
require(flextable)

endpoint_list = c("Normal", "Binary", "Time-to-event")

rowMax = function(x) {

  row_max = rep(0, nrow(x))

  for (i in 1:nrow(x)) row_max[i] = max(x[i, ])

  return(row_max)

}

# Check if the number is an integer
is.wholenumber = function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol

# Arguments:
# parameter: Parameter's value
# n_values: Required number of values
# lower_values: Lower range
# lower_values_sign: Inequality for evaluating the lower range
# upper_values: Upper range
# upper_values_sign: Inequality for evaluating the upper range
# parameter_name: Parameter's name
# component_name: Names of the individual components
# type: Parameter's type (double or integer)
# default_value: Default value
ContinuousErrorCheck = function(parameter, n_values, lower_values, lower_values_sign, upper_values, upper_values_sign, parameter_name, component_name, type = "double", default_value = NA) {

    if (is.null(parameter)) {

        if (!is.na(default_value)) {
            for (i in 1:n_values) {
                parameter[i] = default_value
            }
            return(parameter)
        } else {
            error_message = paste0(parameter_name, " must be specified.") 
            stop(error_message, call. = FALSE)
        }
    } 

    if (!is.na(n_values)) {

      if (length(parameter) != n_values) {
          error_message = paste0(parameter_name, ": ", n_values, " values must be specified.") 
          stop(error_message, call. = FALSE)
      } 

    } else {

      n_values = length(parameter)

    }
    
    if (length(component_name) == 1) {
      if(is.na(component_name[1])) component_name = rep("Each value", n_values)
    }

    for (i in 1:n_values) {

        if (type == "double") {

            if (!is.numeric(parameter[i])) {
                error_message = paste0(parameter_name, ": ", component_name[i], " must be numeric.") 
                stop(error_message, call. = FALSE)            
            }
        }

        if (type == "integer" || type == "int") {

            if (!is.wholenumber(parameter[i])) {
                error_message = paste0(parameter_name, ": ", component_name[i], " must be an integer.") 
                stop(error_message, call. = FALSE)            
            }
        }

        if (length(lower_values) == 1) {

          if (!is.na(lower_values)) {            
              if (lower_values_sign == ">" & parameter[i] <= lower_values) {
                  error_message = paste0(parameter_name, ": Each value must be > ", lower_values, ".") 
                  stop(error_message, call. = FALSE)
              }
              if (lower_values_sign == ">=" & parameter[i] < lower_values) {
                  error_message = paste0(parameter_name, ": Each value must be >= ", lower_values, ".") 
                  stop(error_message, call. = FALSE)
              }
          }


        } else {

          if (!is.na(lower_values[i])) {            
              if (lower_values_sign[i] == ">" & parameter[i] <= lower_values[i]) {
                  error_message = paste0(parameter_name, ": ", component_name[i], " must be > ", lower_values[i], ".") 
                  stop(error_message, call. = FALSE)
              }
              if (lower_values_sign[i] == ">=" & parameter[i] < lower_values[i]) {
                  error_message = paste0(parameter_name, ": ", component_name[i], " must be >= ", lower_values[i], ".") 
                  stop(error_message, call. = FALSE)
              }
          }

        }

        if (length(upper_values) == 1) {

          if (!is.na(upper_values)) {            
              if (upper_values_sign == "<" & parameter[i] >= upper_values) {
                  error_message = paste0(parameter_name, ": Each value must be < ", upper_values, ".") 
                  stop(error_message, call. = FALSE)
              }
              if (upper_values_sign == "<=" & parameter[i] > upper_values) {
                  error_message = paste0(parameter_name, ": Each value must be <= ", upper_values, ".") 
                  stop(error_message, call. = FALSE)
              }
          }

        } else {

          if (!is.na(upper_values[i])) {
              if (upper_values_sign[i] == "<" & parameter[i] >= upper_values[i]) {
                  error_message = paste0(parameter_name, ": ", component_name[i], " must be < ", upper_values[i], ".") 
                  stop(error_message, call. = FALSE)
              }
              if (upper_values_sign[i] == "<=" & parameter[i] > upper_values[i]) {
                  error_message = paste0(parameter_name, ": ", component_name[i], " must be <= ", upper_values[i], ".") 
                  stop(error_message, call. = FALSE)
              }
            }
        }
    }

    return(parameter)
}

# nocov start

ConvertFromCPToHR = function(cp, n_interim, n_final, alpha) {

  w1 = n_interim / n_final
  w2 = 1 - w1
  z_alpha = qnorm(1 - alpha)

  val = qnorm(cp)

  effect_size = (val + z_alpha / sqrt(w2)) / (sqrt((n_final - n_interim) / 4) + sqrt(n_interim / 4) * sqrt(w1 / w2))   
  hr = exp(-effect_size)

  return(hr)

}

# nocov end

SaveReport = function(report, report_title) {

  # Create a docx object
  doc = officer::read_docx(system.file(package = "MedianaDesigner", "template/report_template.docx"))
  dim_doc = officer::docx_dim(doc)

  # Report's title
  doc = officer::set_doc_properties(doc, title = report_title)
  doc = officer::body_add_par(doc, value = report_title, style = "heading 1")

  # Text formatting
  my.text.format = officer::fp_text(font.size = 12, font.family = "Arial")

  # Table formatting
  header.cellProperties = officer::fp_cell(border.left = officer::fp_border(width = 0), border.right = officer::fp_border(width = 0), border.bottom = officer::fp_border(width = 2), border.top = officer::fp_border(width = 2), background.color = "#eeeeee")
  data.cellProperties = officer::fp_cell(border.left = officer::fp_border(width = 0), border.right = officer::fp_border(width = 0), border.bottom = officer::fp_border(width = 0), border.top = officer::fp_border(width = 0))

  header.textProperties = officer::fp_text(font.size = 12, bold = TRUE, font.family = "Arial")
  data.textProperties = officer::fp_text(font.size = 12, font.family = "Arial")

  thick_border = fp_border(color = "black", width = 2)

  leftPar = officer::fp_par(text.align = "left")
  rightPar = officer::fp_par(text.align = "right")
  centerPar = officer::fp_par(text.align = "center")

  # Number of sections in the report (the report's title is not counted)
  n_sections = length(report) 

  # Loop over the sections in the report
  for(section_index in 1:n_sections) {

      # Determine the item's type (text by default)
      type = report[[section_index]]$type

      # Determine the item's label 
      label = report[[section_index]]$label

      # Determine the item's footnote 
      footnote = report[[section_index]]$footnote

      # Determine the item's value 
      value = report[[section_index]]$value

      # Determine column width 
      column_width = report[[section_index]]$column_width

      # Determine the page break status 
      page_break = report[[section_index]]$page_break
      if (is.null(page_break)) page_break = FALSE

      # Determine the figure's location (for figures only)
      filename = report[[section_index]]$filename

      # Determine the figure's dimensions (for figures only)
      dim = report[[section_index]]$dim

      if (!is.null(type)) {

        # Single paragraph 
        if (type == "paragraph") {

            doc = officer::body_add_par(doc, value = label, style = "heading 2")

            doc = officer::body_add_par(doc, value = value, style = "No Spacing")

        }

        # Fully formatted data frame 
        if (type == "table") {

            doc = officer::body_add_par(doc, value = label, style = "heading 2")

            summary_table = flextable::regulartable(data = value)
            summary_table = flextable::style(summary_table, pr_p = leftPar, pr_c = header.cellProperties, pr_t = header.textProperties, part = "header")
            summary_table = flextable::style(summary_table, pr_p = leftPar, pr_c = data.cellProperties, pr_t = data.textProperties, part = "body")

            summary_table = flextable::hline_bottom(summary_table, part = "body", border = thick_border )

            summary_table = flextable::width(summary_table, width = column_width)

            doc = flextable::body_add_flextable(doc, summary_table)

            if (!is.null(footnote)) doc = officer::body_add_par(doc, value = footnote, style = "Normal")

            if (page_break) doc = officer::body_add_break(doc, pos = "after")

        }

        # Enhanced metafile graphics produced by package devEMF 
        if (type == "emf_plot") {

            doc = officer::body_add_par(doc, value = label, style = "heading 2")

            doc = officer::body_add_img(doc, src = filename, width = dim[1], height = dim[2]) 

            if (!is.null(footnote)) doc = officer::body_add_par(doc, value = footnote, style = "Normal")

            if (page_break) doc = officer::body_add_break(doc, pos = "after")

            # Delete the figure
            if (file.exists(filename)) file.remove(filename)   

        }

      }    

  }

  return(doc)          

}
# End SaveReport

CreateTable = function(data_frame, column_names, column_width, title, page_break, footnote = NULL) {

    if (is.null(column_width)) {
      column_width = rep(2, dim(data_frame)[2])
    } 
     
    data_frame = as.data.frame(data_frame)

    data_frame = data.frame(lapply(data_frame, as.character), stringsAsFactors = FALSE)

    colnames(data_frame) = column_names

    item_list = list(label = title, 
                     value = data_frame,
                     column_width = column_width,
                     type = "table",
                     footnote = footnote,
                     page_break = page_break)

    return(item_list)

}
# End of CreateTable       

# Generate function-specific simulation reports 
ReportDoc = function(results) {

  if (class(results) == "ADSSModResults") doc = ADSSModReportDoc(results)
  if (class(results) == "ADTreatSelResults") doc = ADTreatSelReportDoc(results)
  if (class(results) == "ADPopSelResults") doc = ADPopSelReportDoc(results)
  if (class(results) == "FutRuleResults") doc = FutRuleReportDoc(results)
  if (class(results) == "EventPredResults") doc = EventPredReportDoc(results)
  if (class(results) == "ADRandResults") doc = ADRandReportDoc(results)
  if (class(results) == "MultAdjResults") doc = MultAdjReportDoc(results)

  return(doc)  

}
# End of ReportDoc       

# Generate and print function-specific simulation reports 
GenerateReport = function(results, report_filename) {

  # Print the report
  print(ReportDoc(results), target = report_filename)          

}
# End of GenerateReport       

#' Print method for ADSSModResults
#'
#' @param x ADSSModResults object
#' @param ... Arguments passed to or from other methods
#'
#' @export
#' @exportS3Method
print.ADSSModResults = function (x, ...) {
  cat("Use the GenerateReport function to create a detailed simulation report.\n")  # nocov
}

#' Print method for ADTreatSelResults
#'
#' @param x ADTreatSelResults object
#' @param ... Arguments passed to or from other methods
#'
#' @export
#' @exportS3Method
print.ADTreatSelResults = function (x, ...) {
  cat("Use the GenerateReport function to create a detailed simulation report.\n")  # nocov
}

#' Print method for ADPopSelResults
#'
#' @param x ADPopSelResults object
#' @param ... Arguments passed to or from other methods
#'
#' @export
#' @exportS3Method
print.ADPopSelResults = function (x, ...) {
  cat("Use the GenerateReport function to create a detailed simulation report.\n")  # nocov
}

#' Print method for FutRuleResults
#'
#' @param x FutRuleResults object
#' @param ... Arguments passed to or from other methods
#'
#' @export
#' @exportS3Method
print.FutRuleResults = function (x, ...) {
  cat("Use the GenerateReport function to create a detailed simulation report.\n")  # nocov
}

#' Print method for EventPredResults
#'
#' @param x EventPredResults object
#' @param ... Arguments passed to or from other methods
#'
#' @export
#' @exportS3Method
print.EventPredResults = function (x, ...) {
  cat("Use the GenerateReport function to create a detailed simulation report.\n")  # nocov
}

#' Print method for ADRandResults
#'
#' @param x ADRandResults object
#' @param ... Arguments passed to or from other methods
#'
#' @export
#' @exportS3Method
print.ADRandResults = function (x, ...) {
  cat("Use the GenerateReport function to create a detailed simulation report.\n")  # nocov
}

#' Print method for MultAdjResults
#'
#' @param x MultAdjResults object
#' @param ... Arguments passed to or from other methods
#'
#' @export
#' @exportS3Method
print.MultAdjResults = function (x, ...) {
  cat("Use the GenerateReport function to create a detailed simulation report.\n")  # nocov
}
