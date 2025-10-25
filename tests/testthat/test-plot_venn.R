library(testthat)
library(damidBind)

context("Visualization: plot_venn")

# Dummy diff_results object for testing
make_dummy_diff_results_for_venn <- function() {
    analysis_table <- data.frame(
        logFC = c(2.5, -3.0, 0.5, 1.8, -1.2, 0.1, 0.01),
        adj.P.Val = c(0.001, 0.005, 0.1, 0.002, 0.08, 0.5, 0.2), # L1,L2,L4 sig; L5 non-sig; L3,L6 non-sig
        row.names = c("Locus1", "Locus2", "Locus3", "Locus4", "Locus5", "Locus6", "Locus7")
    )

    upCond1 <- analysis_table[c("Locus1", "Locus4"), ] # Sig for Cond1
    upCond2 <- analysis_table[c("Locus2"), ] # Sig for Cond2

    new("DamIDResults",
        upCond1 = upCond1,
        upCond2 = upCond2,
        analysis = analysis_table,
        cond = c(CondA = "Treatment_A", CondB = "Treatment_B"),
        data = list(test_category = "bound")
    )
}

test_that("plot_venn correctly prepares data and calls BioVenn::draw.venn", {
    diff_res <- make_dummy_diff_results_for_venn()

    # Loci: Locus1,2,3,4,5,6,7
    # upCond1 (CondA): Locus1, Locus4
    # upCond2 (CondB): Locus2
    # All unique sig: Locus1, Locus2, Locus4
    # Non-sig: Locus3, Locus5, Locus6, Locus7

    # Cond1_full: union(upCond1, nonsig) = {Locus1, Locus4, Locus3, Locus5, Locus6, Locus7}
    # Cond2_full: union(upCond2, nonsig) = {Locus2, Locus3, Locus5, Locus6, Locus7}

    # Define a mock for BioVenn::draw.venn
    # This mock will store the parameters it was called with.
    captured_args <- new.env()
    mock_draw_venn <- function(list_x, list_y, ..., filename = NULL, output = NULL) {
        captured_args$last_call <- list(
            list_x = list_x,
            list_y = list_y,
            filename = filename,
            output = output,
            other_args = list(...)
        )
        message("BioVenn::draw.venn mocked successfully")
        invisible(NULL)
    }

    local_mocked_bindings(draw.venn = mock_draw_venn, .package = "BioVenn")

    # Define specific set labels for the test
    custom_set_labels <- c("Set A Features", "Set B Features")

    # Call plot_venn
    expect_no_error({
        plot_venn(
            diff_res,
            title = "My Venn Plot",
            subtitle = "Test Subtitle",
            set_labels = custom_set_labels,
            filename = "test_venn.pdf",
            format = "pdf"
        )
    })

    # Verify mock was called
    expect_true(exists("last_call", envir = captured_args))

    # Expected non-significant loci:
    nonsig_expected <- setdiff(rownames(analysisTable(diff_res)), union(rownames(enrichedCond1(diff_res)), rownames(enrichedCond2(diff_res))))
    expect_equal(sort(captured_args$last_call$list_x), sort(union(rownames(enrichedCond1(diff_res)), nonsig_expected)))
    expect_equal(sort(captured_args$last_call$list_y), sort(union(rownames(enrichedCond2(diff_res)), nonsig_expected)))

    # Verify other parameters
    expect_equal(captured_args$last_call$other_args$xtitle, custom_set_labels[1])
    expect_equal(captured_args$last_call$other_args$ytitle, custom_set_labels[2])
    expect_equal(captured_args$last_call$other_args$title, "My Venn Plot")
    expect_equal(captured_args$last_call$other_args$subtitle, "Test Subtitle")
    expect_equal(captured_args$last_call$filename, "test_venn.pdf")
    expect_equal(captured_args$last_call$output, "pdf")
})

test_that("plot_venn handles cases with no significant regions", {
    base_results <- make_dummy_diff_results_for_venn()

    # dummy results object with no significant hits
    diff_res_no_sig <- new("DamIDResults",
        upCond1 = analysisTable(base_results)[FALSE, ],
        upCond2 = analysisTable(base_results)[FALSE, ],
        analysis = analysisTable(base_results),
        cond = conditionNames(base_results),
        data = inputData(base_results)
    )

    mock_draw_venn_args <- NULL # Reset mock args
    mock_draw_venn <- function(list_x, list_y, ...) {
        mock_draw_venn_args <<- list(list_x = list_x, list_y = list_y)
        invisible(NULL)
    }
    local_mocked_bindings(draw.venn = mock_draw_venn, .package = "BioVenn")

    res <- evaluate_promise(plot_venn(diff_res_no_sig))
    expect_length(res$messages, 2)
    expect_match(res$messages[1], "Note: No loci were differentially enriched in condition 1", fixed = TRUE)
    expect_match(res$messages[2], "Note: No loci were differentially enriched in condition 2", fixed = TRUE)

    expect_true(!is.null(mock_draw_venn_args))
    # All loci should be in both sets, as all are non-significant
    expect_equal(sort(mock_draw_venn_args$list_x), sort(rownames(analysisTable(diff_res_no_sig))))
    expect_equal(sort(mock_draw_venn_args$list_y), sort(rownames(analysisTable(diff_res_no_sig))))
})


test_that("plot_venn messages if one condition has no loci", {
    base_results <- make_dummy_diff_results_for_venn()

    # dummy results object with no significant hits
    diff_res_partial_sig <- new("DamIDResults",
        upCond1 = analysisTable(base_results)[FALSE, ],
        upCond2 = enrichedCond2(base_results),
        analysis = analysisTable(base_results),
        cond = conditionNames(base_results),
        data = inputData(base_results)
    )

    mock_draw_venn <- function(list_x, list_y, ...) invisible(NULL)
    local_mocked_bindings(draw.venn = mock_draw_venn, .package = "BioVenn")

    expect_message(plot_venn(diff_res_partial_sig), "Note: No loci were differentially enriched in condition 1")
})

test_that("plot_venn uses default set_labels if not provided", {
    diff_res <- make_dummy_diff_results_for_venn()
    captured_args_env <- new.env()
    mock_draw_venn <- function(list_x, list_y, xtitle, ytitle, ...) {
        captured_args_env$venn_args <- list(xtitle = xtitle, ytitle = ytitle)
        invisible(NULL)
    }
    local_mocked_bindings(draw.venn = mock_draw_venn, .package = "BioVenn")

    expect_no_error(plot_venn(diff_res, set_labels = NULL))
    expect_true(!is.null(captured_args_env$venn_args))
    expect_equal(captured_args_env$venn_args$xtitle, names(conditionNames(diff_res))[1])
    expect_equal(captured_args_env$venn_args$ytitle, names(conditionNames(diff_res))[2])
})
