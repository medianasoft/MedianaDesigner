test_that("CreateTable works", {
  
  check = function(column_width) {

    column_names = c("Trial arm", "Sample size")
    trial_arms = c("Control", "Treatment")
    sample_size = 10
    data_frame = data.frame(trial_arms, sample_size)
    title = paste0("Table 1. Number of enrolled patients")

    column_width = column_width
    table = CreateTable(data_frame, column_names, column_width, title, FALSE)

    expect_equal(table$label, title)

    t = table$value
    expect_equal(typeof(t), "list")
    expect_equal(length(t), 2)
    expect_equal(colnames(t), column_names)
    expect_equal(rownames(t), c("1","2"))

    if (is.null(column_width)) {
      expect_equal(table$column_width, c(2,2))
    } else {
      expect_equal(table$column_width, column_width)
    }

    expect_equal(table$type, "table")
    expect_equal(table$footnote, NULL)
    expect_equal(table$page_break, FALSE)

  }

  check(c(5, 1.5))
  check(NULL)
})
