# Universal test-function

testParameterErrors = function(params, paramName,
    checkMissing = FALSE, checkWrong = NA, checkSize = FALSE, checkMin = NA, checkMax = NA) {

    func = ClustRand

    isParam = function(param) {
        return( length(param)==1 && !is.na(param) || length(param)>1 )
    }

    # Missing
    if (checkMissing) {
        testParams = params
        testParams[paramName] <- NULL
        expect_error(func(testParams),
            info = paste0("Checking for missing ", paramName))
    }
    # Wrong
    if (isParam(checkWrong)) {
        testParams = params
        testParams[[paramName]] <- checkWrong
        expect_error(func(testParams),
            info = paste0("Checking for wrong ", paramName))
    }
    # Check size
    if (checkSize) {
        testParams = params
        testParams[[paramName]] <- append(testParams[[paramName]], testParams[[paramName]][1])
        expect_error(func(testParams),
            info = paste0("Checking for value size ", paramName))
    }
    # Check below min value
    if (isParam(checkMin)) {
        testParams = params
        testParams[[paramName]] <- checkMin
        expect_error(func(testParams),
            info = paste0("Checking for min value ", paramName))
    }
    # Check under max value
    if (isParam(checkMax)) {
        testParams = params
        testParams[[paramName]] <- checkMax
        expect_error(func(testParams),
            info = paste0("Checking for max value ", paramName))
    }
}
